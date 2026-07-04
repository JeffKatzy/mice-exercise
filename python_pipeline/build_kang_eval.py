"""
Full pipeline for Kang et al. 2018 (GSE96583) eval set.
Downloads batch2 count matrix, assembles h5ad, runs DESeq2 pseudobulk
(stim vs ctrl per cell type), writes sig-only CSV.

Run once; subsequent runs skip already-completed steps.
"""

import urllib.request
import gzip
import shutil
import tarfile
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.io
import anndata as ad

OUT_DIR   = Path(__file__).parent / "dataset" / "kang_pbmc"
OUTPUTS   = Path(__file__).parent / "outputs"
BASE_FTP  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── Step 1: Download ──────────────────────────────────────────────────────────

def download(url: str, dest: Path):
    if dest.exists():
        print(f"  skip (exists): {dest.name}")
        return
    print(f"  downloading {dest.name}...")
    urllib.request.urlretrieve(url, dest)
    print(f"  done ({dest.stat().st_size // 1024} KB)")


def decompress_gz(src: Path, dest: Path):
    if dest.exists():
        return
    with gzip.open(src, 'rb') as f_in, open(dest, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def step1_download():
    print("\n── Step 1: Download ──")
    raw_tar = OUT_DIR / "GSE96583_RAW.tar"
    meta_gz = OUT_DIR / "GSE96583_batch2.total.tsne.df.tsv.gz"
    genes_gz = OUT_DIR / "GSE96583_batch2.genes.tsv.gz"

    download(f"{BASE_FTP}/GSE96583_RAW.tar", raw_tar)
    download(f"{BASE_FTP}/GSE96583_batch2.total.tsne.df.tsv.gz", meta_gz)
    download(f"{BASE_FTP}/GSE96583_batch2.genes.tsv.gz", genes_gz)

    # Extract only batch2 files from RAW.tar (GSM2560248=ctrl, GSM2560249=stim)
    batch2_dir = OUT_DIR / "batch2"
    batch2_dir.mkdir(exist_ok=True)
    marker = batch2_dir / ".extracted"
    if not marker.exists():
        print("  extracting batch2 from RAW.tar...")
        with tarfile.open(raw_tar) as tar:
            for member in tar.getmembers():
                if "GSM2560248" in member.name or "GSM2560249" in member.name:
                    member.name = Path(member.name).name
                    tar.extract(member, batch2_dir, filter='data')
        marker.touch()
        print(f"  extracted: {[f.name for f in batch2_dir.iterdir() if not f.name.startswith('.')]}")

    return batch2_dir, meta_gz, genes_gz


# ── Step 2: Assemble h5ad ────────────────────────────────────────────────────

def step2_assemble(batch2_dir: Path, meta_gz: Path, genes_gz: Path) -> Path:
    h5ad_path = OUT_DIR / "kang_pbmc.h5ad"
    if h5ad_path.exists():
        print(f"\n── Step 2: Assemble — skip (exists: {h5ad_path.name}) ──")
        return h5ad_path

    print("\n── Step 2: Assemble h5ad ──")

    # Load metadata
    meta_tsv = OUT_DIR / "GSE96583_batch2.total.tsne.df.tsv"
    decompress_gz(meta_gz, meta_tsv)
    meta = pd.read_csv(meta_tsv, sep="\t", index_col=0)

    # Load shared genes list (one file covers both samples)
    genes_tsv = OUT_DIR / "GSE96583_batch2.genes.tsv"
    decompress_gz(genes_gz, genes_tsv)
    genes = pd.read_csv(genes_tsv, header=None)[0].tolist()

    # GSM2560248 = ctrl (2.1), GSM2560249 = stim (2.2)
    sample_map = {
        'GSM2560248': 'ctrl',
        'GSM2560249': 'stim',
    }

    mtx_files = sorted(batch2_dir.glob("*.mtx.gz"))
    barcodes_files = list(batch2_dir.glob("*barcodes*"))

    print(f"  matrix files: {[f.name for f in mtx_files]}")
    print(f"  barcodes files: {[f.name for f in barcodes_files]}")

    adatas = []
    for mtx_file in mtx_files:
        gsm = mtx_file.name.split("_")[0]
        bc_file = next((f for f in barcodes_files if gsm in f.name), None)
        if bc_file is None:
            print(f"  WARNING: no barcodes for {gsm}")
            continue

        mat = scipy.io.mmread(gzip.open(mtx_file, 'rb'))
        mat = sp.csr_matrix(mat.T)

        bc = pd.read_csv(bc_file, header=None, compression='gzip')[0].tolist()
        adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=genes))
        adata.obs['sample'] = gsm
        adata.obs['sample_stim'] = sample_map.get(gsm, 'unknown')
        adatas.append(adata)
        print(f"  {gsm} ({sample_map.get(gsm,'?')}): {adata.shape}")

    adata = ad.concat(adatas, join='outer')
    adata.obs_names_make_unique()

    # Join metadata — match on barcode index
    adata.obs = adata.obs.join(meta[['stim', 'cell', 'ind', 'multiplets']], how='left')

    # Filter singlets only
    adata = adata[adata.obs['multiplets'] == 'singlet'].copy()
    adata = adata[adata.obs['stim'].notna() & adata.obs['cell'].notna()].copy()

    print(f"  final shape: {adata.shape}")
    print(f"  stim: {adata.obs['stim'].value_counts().to_dict()}")
    print(f"  cell types: {adata.obs['cell'].value_counts().to_dict()}")

    adata.write_h5ad(h5ad_path)
    print(f"  saved: {h5ad_path}")
    return h5ad_path


# ── Step 3: Pseudobulk DESeq2 ────────────────────────────────────────────────

def step3_degs(h5ad_path: Path) -> Path:
    out_csv = OUTPUTS / "degs_kang_pbmc_deseq2_pseudobulk.csv"
    sig_csv = OUTPUTS / "degs_kang_pbmc_deseq2_sig_only.csv"
    if sig_csv.exists():
        print(f"\n── Step 3: DEGs — skip (exists: {sig_csv.name}) ──")
        return sig_csv

    print("\n── Step 3: Pseudobulk DESeq2 ──")
    import scanpy as sc
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    adata = sc.read_h5ad(h5ad_path)
    adata.var_names_make_unique()

    results = []

    for cell_type in adata.obs['cell'].unique():
        sub = adata[adata.obs['cell'] == cell_type].copy()

        # Pseudobulk: sum counts per donor per condition
        donors = sub.obs['ind'].astype(str)
        conditions = sub.obs['stim']
        pseudo_rows = []
        labels = []

        for ind in sub.obs['ind'].unique():
            for stim in ['ctrl', 'stim']:
                mask = (sub.obs['ind'] == ind) & (sub.obs['stim'] == stim)
                if mask.sum() < 10:
                    continue
                counts = np.asarray(sub[mask].X.sum(axis=0)).ravel()
                pseudo_rows.append(counts)
                labels.append({'ind': str(ind), 'stim': stim})

        if len(pseudo_rows) < 4:
            print(f"  {cell_type}: too few pseudobulk samples ({len(pseudo_rows)}), skipping")
            continue

        pb_df = pd.DataFrame(
            np.array(pseudo_rows),
            columns=sub.var_names,
            index=[f"{r['ind']}_{r['stim']}" for r in labels]
        ).astype(int)
        meta_df = pd.DataFrame(labels, index=pb_df.index)

        # Filter low-count genes
        keep = (pb_df > 0).sum(axis=0) >= 3
        pb_df = pb_df.loc[:, keep]

        n_ctrl = (meta_df['stim'] == 'ctrl').sum()
        n_stim = (meta_df['stim'] == 'stim').sum()

        if n_ctrl < 2 or n_stim < 2:
            print(f"  {cell_type}: n_ctrl={n_ctrl} n_stim={n_stim}, skipping")
            continue

        print(f"  {cell_type}: {pb_df.shape[1]} genes, {n_ctrl} ctrl / {n_stim} stim samples")

        try:
            dds = DeseqDataSet(
                counts=pb_df,
                metadata=meta_df,
                design_factors="stim",
                ref_level=["stim", "ctrl"],
                refit_cooks=True,
                quiet=True,
            )
            dds.deseq2()
            stat_res = DeseqStats(dds, contrast=["stim", "stim", "ctrl"], quiet=True)
            stat_res.summary()
            res = stat_res.results_df.copy()
            res['gene'] = res.index
            res['cell_state'] = cell_type
            res['tissue'] = 'PBMC'
            res['contrast'] = 'stim_vs_ctrl'
            res['contrast_label'] = 'IFNb_stim_vs_ctrl'
            res['dataset'] = 'kang_pbmc'
            res['n_ctrl'] = n_ctrl
            res['n_stim'] = n_stim
            res['is_sig'] = (res['padj'] < 0.05) & (res['padj'].notna())
            results.append(res)
        except Exception as e:
            print(f"  {cell_type}: DESeq2 error — {e}")

    if not results:
        print("  No results computed.")
        return sig_csv

    full = pd.concat(results, ignore_index=True)
    full.to_csv(out_csv, index=False)
    full[full['is_sig']].to_csv(sig_csv, index=False)

    print(f"\nSig DEGs per cell type:")
    print(full[full['is_sig']].groupby('cell_state').size().sort_values(ascending=False).to_string())
    print(f"\nSaved: {out_csv.name} ({len(full)} rows)")
    print(f"Saved: {sig_csv.name} ({full['is_sig'].sum()} sig rows)")
    return sig_csv


if __name__ == "__main__":
    batch2_dir, meta_gz, genes_gz = step1_download()
    h5ad_path = step2_assemble(batch2_dir, meta_gz, genes_gz)
    step3_degs(h5ad_path)
