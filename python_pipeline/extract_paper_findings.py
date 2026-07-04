"""
extract_paper_findings.py

One-time extraction: reads a paper PDF and outputs a structured JSON file of
quantitative findings. This JSON becomes the ground truth for blind test scoring.

Usage:
    python extract_paper_findings.py <paper.pdf> <output.json> [--geo GSE123456]

The JSON is human-reviewable and editable before use as ground truth.
Run once per benchmark paper; scoring against it is then fully automated.
"""

import sys
import os
import json
import base64
import argparse
from pathlib import Path

# Load .env from project root or current dir if present
for _env_path in [Path(__file__).parent / '.env',
                  Path(__file__).parent.parent / '.env',
                  Path('.env')]:
    if _env_path.exists():
        for line in _env_path.read_text().splitlines():
            if line.strip() and not line.startswith('#') and '=' in line:
                k, v = line.split('=', 1)
                os.environ.setdefault(k.strip(), v.strip())
        break

if not os.environ.get('ANTHROPIC_API_KEY'):
    print("ERROR: ANTHROPIC_API_KEY not set.")
    print("Either:")
    print("  export ANTHROPIC_API_KEY=sk-ant-...")
    print("  or create a .env file with: ANTHROPIC_API_KEY=sk-ant-...")
    sys.exit(1)

import anthropic

EXTRACTION_PROMPT = """You are extracting quantitative findings from a scientific paper
for use as ground truth in an automated evaluation benchmark.

The paper is a single-cell RNA-seq study. Extract ALL quantitative findings that a
re-analysis pipeline should be able to reproduce. Focus on:

1. Differentially expressed genes — specific genes that are significantly up or down
   in specific cell types under specific conditions
2. Cell-type-level observations — which cell types respond most to each condition
3. Enrichment results — pathway, gene set, or TF enrichments with statistics
4. Negative findings — cell types or genes explicitly shown NOT to respond

For each finding, extract:
- gene: gene symbol (e.g. "Fndc5", "Nr1d1") — or null if it's a cell-type-level finding
- cell_type: the specific cell type or state (e.g. "Areg", "WAT_IPC", "committed preadipocyte")
- tissue: the tissue (e.g. "vWAT", "scWAT", "skeletal muscle")
- direction: "up", "down", "no_change", or null if not applicable
- condition_a: the comparison group (e.g. "TC", "exercise", "lean")
- condition_b: the reference group (e.g. "SC", "sedentary", "obese")
- statistic: the reported statistic (e.g. "log2FC=1.5, padj=0.008") — quote exactly
- finding_type: "DEG" | "cell_type_response" | "enrichment" | "negative" | "other"
- confidence: "high" (explicitly stated with statistics) | "medium" (mentioned but no stats) | "low" (implied)
- quote: a short direct quote from the paper supporting this finding (1 sentence max)
- notes: any caveats, e.g. "only in rescue contrast" or "trending, not significant"

Return ONLY valid JSON in this format:
{
  "paper_title": "...",
  "paper_doi": "...",
  "geo_accession": "...",
  "species": "mouse" | "human" | "both",
  "tissues": ["vWAT", "scWAT", ...],
  "conditions": ["SC", "TC", "SH", "TH"],
  "cell_types_annotated": ["Areg", "WAT_IPC", ...],
  "findings": [
    {
      "id": 1,
      "gene": "Nr1d2",
      "cell_type": "WAT_IPC",
      "tissue": "vWAT",
      "direction": "up",
      "condition_a": "TC",
      "condition_b": "SC",
      "statistic": "log2FC=1.3, padj=0.008",
      "finding_type": "DEG",
      "confidence": "high",
      "quote": "Nr1d2 was significantly upregulated in WAT_IPC cells...",
      "notes": ""
    },
    ...
  ],
  "key_claims": [
    "MSCs are the central responders to both obesity and exercise",
    "Exercise reactivates circadian clock genes in fat MSCs"
  ],
  "extraction_notes": "Any caveats about extraction quality or ambiguous cases"
}

Be exhaustive — include findings even if they seem minor. A pipeline that finds a
'minor' finding the paper reported is a true positive; missing it is a false negative.
Include negative findings (e.g. 'gene X not significant in cell type Y') as these are
important true negatives for the eval harness."""


def read_pdf_as_base64(pdf_path: str) -> str:
    with open(pdf_path, 'rb') as f:
        return base64.standard_b64encode(f.read()).decode('utf-8')


def extract_findings(pdf_path: str, geo_accession: str = None) -> dict:
    client = anthropic.Anthropic()

    print(f"Reading {pdf_path}...")
    pdf_b64 = read_pdf_as_base64(pdf_path)

    geo_note = f"\nGEO accession for this dataset: {geo_accession}" if geo_accession else ""

    print("Sending to Claude claude-opus-4-7 for extraction...")
    response = client.messages.create(
        model="claude-opus-4-7",
        max_tokens=8192,
        messages=[
            {
                "role": "user",
                "content": [
                    {
                        "type": "document",
                        "source": {
                            "type": "base64",
                            "media_type": "application/pdf",
                            "data": pdf_b64,
                        },
                    },
                    {
                        "type": "text",
                        "text": EXTRACTION_PROMPT + geo_note,
                    }
                ],
            }
        ],
    )

    raw = response.content[0].text.strip()

    # Strip markdown code fences if present
    if raw.startswith("```"):
        lines = raw.split('\n')
        raw = '\n'.join(lines[1:-1] if lines[-1] == '```' else lines[1:])

    try:
        findings = json.loads(raw)
    except json.JSONDecodeError as e:
        print(f"JSON parse error: {e}")
        print("Raw response saved to extraction_raw.txt")
        Path("extraction_raw.txt").write_text(raw)
        raise

    findings['extraction_model'] = "claude-opus-4-7"
    findings['source_pdf'] = str(pdf_path)
    return findings


def score_pipeline_findings(pipeline_findings: list[dict], ground_truth: dict) -> dict:
    """
    Compare pipeline output against ground truth JSON.

    pipeline_findings: list of {gene, cell_type, tissue, direction, contrast}
    ground_truth: output of extract_findings()

    Returns: {tp, fp, fn, precision, recall, details}
    """
    gt_findings = ground_truth['findings']
    high_confidence = [f for f in gt_findings if f['confidence'] == 'high']

    tp, fp, fn_list = 0, 0, []

    def matches(pipeline_f, gt_f):
        gene_match = (
            pipeline_f.get('gene', '').lower() == (gt_f.get('gene') or '').lower()
        )
        cell_match = (
            pipeline_f.get('cell_type', '').lower() in gt_f.get('cell_type', '').lower()
            or gt_f.get('cell_type', '').lower() in pipeline_f.get('cell_type', '').lower()
        )
        dir_match = pipeline_f.get('direction') == gt_f.get('direction')
        return gene_match and cell_match and dir_match

    matched_gt = set()
    for pf in pipeline_findings:
        found = False
        for i, gf in enumerate(high_confidence):
            if matches(pf, gf):
                tp += 1
                matched_gt.add(i)
                found = True
                break
        if not found:
            fp += 1

    for i, gf in enumerate(high_confidence):
        if i not in matched_gt:
            fn_list.append(gf)
            fn += 1 if False else 0  # counted below

    fn = len(high_confidence) - len(matched_gt)
    precision = tp / (tp + fp) if (tp + fp) > 0 else None
    recall = tp / (tp + fn) if (tp + fn) > 0 else None

    return {
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'precision': round(precision, 3) if precision is not None else None,
        'recall': round(recall, 3) if recall is not None else None,
        'missed_findings': fn_list,
        'n_ground_truth_high_confidence': len(high_confidence),
        'n_pipeline_findings': len(pipeline_findings),
    }


def main():
    parser = argparse.ArgumentParser(description='Extract paper findings to JSON benchmark')
    parser.add_argument('pdf', help='Path to paper PDF')
    parser.add_argument('output', help='Output JSON path')
    parser.add_argument('--geo', default=None, help='GEO accession number (e.g. GSE183288)')
    args = parser.parse_args()

    findings = extract_findings(args.pdf, args.geo)

    out_path = Path(args.output)
    out_path.write_text(json.dumps(findings, indent=2))

    print(f"\nExtracted {len(findings['findings'])} findings")
    print(f"  High confidence: {sum(1 for f in findings['findings'] if f['confidence'] == 'high')}")
    print(f"  Medium confidence: {sum(1 for f in findings['findings'] if f['confidence'] == 'medium')}")
    print(f"  Negative findings: {sum(1 for f in findings['findings'] if f['finding_type'] == 'negative')}")
    print(f"\nKey claims ({len(findings.get('key_claims', []))}):")
    for c in findings.get('key_claims', []):
        print(f"  - {c}")
    print(f"\nSaved to: {out_path}")
    print("\nReview and edit the JSON before using as ground truth.")
    print("Especially check: are negative findings complete? Are cell type names consistent?")


if __name__ == '__main__':
    main()
