"""
Build eval set from Jerby-Arnon et al. 2018 Cell (GSE115978).
Melanoma scRNA-seq, responder vs non-responder to anti-PD1 immunotherapy.
Lab: Regev/Hacohen (Broad/MIT). GitHub: livnatje/ImmuneResistance.

Steps:
  1. Download count matrix + cell annotations via GEOparse
  2. Assemble h5ad with cell type + response metadata
  3. Run pseudobulk DESeq2 per cell type (responder vs non-responder)
  4. Write sig-only CSV for eval_degs.py
"""

import GEOparse
import urllib.request
import gzip
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import scipy.sparse as sp
import anndata as ad

GEO_ID   = "GSE115978"
OUT_DIR  = Path(__file__).parent / "dataset" / "jerby_melanoma"
OUTPUTS  = Path(__file__).parent / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BASE_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl"


def decompress_gz(src: Path, dest: Path):
    if dest.exists():
        return
    with gzip.open(src, 'rb') as f_in, open(dest, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def download(url: str, dest: Path):
    if dest.exists():
        print(f"  skip (exists): {dest.name}")
        return
    print(f"  downloading {dest.name}...")
    urllib.request.urlretrieve(url, dest)
    print(f"  done ({dest.stat().st_size // 1024} KB)")


# ── Step 1: Download ──────────────────────────────────────────────────────────

def step1_download():
    print("\n── Step 1: Download ──")

    counts_gz = OUT_DIR / "GSE115978_counts.csv.gz"
    annot_gz  = OUT_DIR / "GSE115978_cell.annotations.csv.gz"

    download(f"{BASE_FTP}/GSE115978_counts.csv.gz", counts_gz)
    download(f"{BASE_FTP}/GSE115978_cell.annotations.csv.gz", annot_gz)

    # Use GEOparse to get per-sample metadata (cell type, treatment, cohort)
    soft_cache = OUT_DIR / f"{GEO_ID}_family.soft.gz"
    if soft_cache.exists():
        print(f"  skip (exists): {soft_cache.name}")
        gse = GEOparse.get_GEO(GEO_ID, destdir=str(OUT_DIR), silent=True)
    else:
        print("  downloading SOFT metadata via GEOparse...")
        gse = GEOparse.get_GEO(GEO_ID, destdir=str(OUT_DIR), silent=True)

    return counts_gz, annot_gz, gse


# ── Step 2: Assemble h5ad ────────────────────────────────────────────────────

def step2_assemble(counts_gz: Path, annot_gz: Path, gse) -> Path:
    h5ad_path = OUT_DIR / "jerby_melanoma.h5ad"
    if h5ad_path.exists():
        print(f"\n── Step 2: Assemble — skip (exists: {h5ad_path.name}) ──")
        return h5ad_path

    print("\n── Step 2: Assemble h5ad ──")

    # Load cell annotations (cell barcode → cell type, patient, response)
    annot_csv = OUT_DIR / "GSE115978_cell.annotations.csv"
    decompress_gz(annot_gz, annot_csv)
    annot = pd.read_csv(annot_csv, index_col=0)
    print(f"  annotations shape: {annot.shape}")
    print(f"  annotation columns: {annot.columns.tolist()}")
    print(f"  value counts:\n{annot.apply(lambda c: c.nunique()).to_string()}")

    # Load count matrix (genes x cells)
    print("  loading count matrix (this may take a minute)...")
    counts_csv = OUT_DIR / "GSE115978_counts.csv"
    decompress_gz(counts_gz, counts_csv)
    counts = pd.read_csv(counts_csv, index_col=0)
    print(f"  counts shape: {counts.shape} (genes x cells)")

    # Align cells
    common_cells = annot.index.intersection(counts.columns)
    print(f"  common cells: {len(common_cells)}")
    counts = counts[common_cells]
    annot  = annot.loc[common_cells]

    # Build AnnData (cells x genes)
    mat = sp.csr_matrix(counts.values.T)
    adata = ad.AnnData(
        X=mat,
        obs=annot,
        var=pd.DataFrame(index=counts.index),
    )
    print(f"  AnnData shape: {adata.shape}")
    print(f"\n  obs columns: {adata.obs.columns.tolist()}")
    for col in adata.obs.columns:
        n = adata.obs[col].nunique()
        if n <= 20:
            print(f"    {col}: {adata.obs[col].value_counts().to_dict()}")

    adata.write_h5ad(h5ad_path)
    print(f"\n  saved: {h5ad_path}")
    return h5ad_path


# ── Step 3: Pseudobulk DESeq2 ────────────────────────────────────────────────

def step3_degs(h5ad_path: Path) -> Path:
    sig_csv = OUTPUTS / "degs_jerby_melanoma_deseq2_sig_only.csv"
    if sig_csv.exists():
        print(f"\n── Step 3: DEGs — skip (exists: {sig_csv.name}) ──")
        return sig_csv

    print("\n── Step 3: Pseudobulk DESeq2 ──")
    import scanpy as sc
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    adata = sc.read_h5ad(h5ad_path)
    adata.var_names_make_unique()

    print(f"  obs columns: {adata.obs.columns.tolist()}")

    # Identify the response and cell type columns from what's available
    # Will print and exit if columns are unexpected so we can adjust
    print(f"  obs sample:\n{adata.obs.head(3).to_string()}")

    # Columns confirmed from annotations file:
    # cell.types, treatment.group (treatment.naive vs post.treatment), samples (patient)
    cell_type_col = "cell.types"
    response_col  = "treatment.group"
    patient_col   = "samples"

    print(f"\n  using: cell_type={cell_type_col}, response={response_col}, patient={patient_col}")

    # Keep only cells with response annotation
    adata = adata[adata.obs[response_col].notna()].copy()
    print(f"  response values: {adata.obs[response_col].value_counts().to_dict()}")
    print(f"  cell types: {adata.obs[cell_type_col].value_counts().to_dict()}")

    results = []
    for cell_type in adata.obs[cell_type_col].unique():
        sub = adata[adata.obs[cell_type_col] == cell_type].copy()

        # Pseudobulk per patient per response group
        pseudo_rows, labels = [], []
        for patient in sub.obs[patient_col].unique():
            for resp in sub.obs[response_col].unique():
                mask = (sub.obs[patient_col] == patient) & (sub.obs[response_col] == resp)
                if mask.sum() < 5:
                    continue
                counts_arr = np.asarray(sub[mask].X.sum(axis=0)).ravel()
                pseudo_rows.append(counts_arr)
                labels.append({'patient': str(patient), 'response': str(resp)})

        if len(pseudo_rows) < 4:
            print(f"  {cell_type}: too few pseudobulk samples ({len(pseudo_rows)}), skipping")
            continue

        pb_df   = pd.DataFrame(np.array(pseudo_rows).astype(int),
                               columns=sub.var_names,
                               index=[f"{r['patient']}_{r['response']}" for r in labels])
        meta_df = pd.DataFrame(labels, index=pb_df.index)

        keep = (pb_df > 0).sum(axis=0) >= 3
        pb_df = pb_df.loc[:, keep]

        resp_vals = meta_df['response'].unique()
        if len(resp_vals) < 2:
            print(f"  {cell_type}: only one response group, skipping")
            continue

        # Determine reference (non-responder) level
        ref = next((v for v in resp_vals if 'NR' in v or 'non' in v.lower() or 'resist' in v.lower()), resp_vals[0])
        test = next((v for v in resp_vals if v != ref), resp_vals[1])

        n_ref  = (meta_df['response'] == ref).sum()
        n_test = (meta_df['response'] == test).sum()
        print(f"  {cell_type}: {pb_df.shape[1]} genes, {n_ref} {ref} / {n_test} {test}")

        if n_ref < 2 or n_test < 2:
            print(f"    insufficient replicates, skipping")
            continue

        try:
            dds = DeseqDataSet(
                counts=pb_df,
                metadata=meta_df,
                design="~ response",
                quiet=True,
            )
            dds.deseq2()
            stat_res = DeseqStats(dds, contrast=["response", test, ref], quiet=True)
            stat_res.summary()
            res = stat_res.results_df.copy()
            res['gene']           = res.index
            res['cell_state']     = cell_type
            res['tissue']         = 'tumor'
            res['contrast']       = f'{test}_vs_{ref}'
            res['contrast_label'] = 'responder_vs_nonresponder'
            res['dataset']        = 'jerby_melanoma'
            res['is_sig']         = (res['padj'] < 0.05) & res['padj'].notna()
            results.append(res)
        except Exception as e:
            print(f"  {cell_type}: DESeq2 error — {e}")

    if not results:
        print("  No results computed.")
        return sig_csv

    full = pd.concat(results, ignore_index=True)
    full.to_csv(OUTPUTS / "degs_jerby_melanoma_deseq2_pseudobulk.csv", index=False)
    full[full['is_sig']].to_csv(sig_csv, index=False)

    print(f"\nSig DEGs per cell type:")
    print(full[full['is_sig']].groupby('cell_state').size().sort_values(ascending=False).to_string())
    print(f"\nSaved: {sig_csv.name} ({full['is_sig'].sum()} sig rows)")
    return sig_csv


if __name__ == "__main__":
    counts_gz, annot_gz, gse = step1_download()
    h5ad_path = step2_assemble(counts_gz, annot_gz, gse)
    step3_degs(h5ad_path)
