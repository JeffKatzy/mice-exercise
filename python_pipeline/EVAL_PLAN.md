# Agentic Pipeline Eval Plan

## Goal
Build a ground-truth DEG eval set across multiple independent datasets so the agentic
pipeline (NB16) can be scored: given a dataset + condition contrast, does it recover
the known significant genes per cell type?

## Current eval set (`eval_degs.py`)

| Dataset | Tissue | Contrast | Cell types | Source |
|---------|--------|----------|-----------|--------|
| yang_exercise | vWAT (mouse) | exercise_TC_vs_SC | 6 | DESeq2 computed (this project) |
| yang_exercise | vWAT (mouse) | obesity_SH_vs_SC | 6 | DESeq2 computed (this project) |
| yang_exercise | vWAT (mouse) | rescue_TH_vs_SH | 7 | DESeq2 computed (this project) |
| mathys_ad | PFC (human) | AD_vs_control | 6 | paper-reported DEGs (reanalysis_Mathys_2019/) |

---

## Datasets to add

### 1. Kang et al. 2018 — IFN-β stimulated PBMCs
- **GEO:** GSE96583 ✓ open, count matrix available
- **Lab:** Satija lab (Seurat paper, top scRNA-seq methods lab)
- **Biology:** 8 donors, PBMCs, stimulated with IFN-β vs control
- **Cell types:** ~8 immune cell types (CD4 T, CD8 T, NK, B, Monocyte, DC, etc.)
- **Why useful:** Canonical benchmark dataset; DEGs well-established; stimulation
  vs control is a clean contrast with large effect sizes
- **GitHub:** satijalab/seurat tutorials reference this extensively
- **Status:** TODO

### 2. Wilk et al. 2020 — COVID-19 PBMCs
- **GEO:** GSE149689 ✓ open, barcodes + features + matrix all available
- **Lab:** Bhatt/Blish lab (Stanford, Nature Medicine)
- **Biology:** Severe COVID vs mild vs healthy, immune cells
- **Cell types:** PBMC immune subtypes
- **Why useful:** Well-cited disease vs healthy contrast; different from stimulation (Kang)
- **Status:** TODO

### 3. Nathan et al. 2021 — RA synovial tissue  
- **GEO:** GSE159117 ✓ open (RAW.tar — need to check cell type annotation availability)
- **Lab:** Raychaudhuri lab (Harvard/Broad, top immunogenomics lab)
- **Biology:** Rheumatoid arthritis PBMCs, inflamed vs healthy
- **Cell types:** T cells, B cells, monocytes, fibroblasts
- **Why useful:** Stromal + immune mix; fibroblast DEGs especially interesting
  (analogous to Yang vWAT fibroblasts)
- **Status:** TODO — need to verify cell type annotations in supplementary

### 4. Smillie et al. 2019 — Ulcerative colitis colon (Regev lab)
- **GEO:** Not confirmed accessible — GSE125527 is a UC dataset but may not be
  the Smillie Cell paper
- **Lab:** Regev lab (Broad Institute)
- **Status:** BLOCKED — need to confirm correct accession

---

## Pipeline per new dataset

```
1. download_<dataset>.py  — fetch count matrix + metadata from GEO
2. assemble_<dataset>.py  — build h5ad (same pattern as assemble_h5ad.py)
3. run_degs_<dataset>.py  — DESeq2 pseudobulk per cell type, case vs control
4. eval_degs.py           — add _load_<dataset>() and register in load_eval_degs()
```

Steps 1-2 can reuse a shared `geo_download.py` utility (wget from FTP + decompress).
Steps 3 uses the same `run_pseudobulk_de.py` pattern already in this project.

---

## eval_degs.py EvalKey structure

```python
EvalKey(dataset, tissue, contrast, cell_state) -> list[str]
```

Adding `dataset` as first dimension means Yang and all new datasets coexist in one dict.
Filter by dataset for dataset-specific scoring.

---

## Scoring (to implement)

Given pipeline output (predicted DEG list per cell type) vs ground truth:

```
precision = |predicted ∩ truth| / |predicted|
recall    = |predicted ∩ truth| / |truth|
F1        = harmonic mean
```

Report per (dataset, contrast, cell_state) and aggregate.
False positives are the dangerous failure mode — prioritize precision.
