"""
Download and explore Kang et al. 2018 PBMC IFN-beta dataset (GSE96583).
Step 1: download metadata/annotation files only to check schema before pulling full matrix.
"""

import urllib.request
import gzip
import shutil
from pathlib import Path
import pandas as pd

OUT_DIR = Path(__file__).parent / "dataset" / "kang_pbmc"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl"

FILES = {
    "batch1_tsne": f"{BASE}/GSE96583_batch1.total.tsne.df.tsv.gz",
    "batch2_tsne": f"{BASE}/GSE96583_batch2.total.tsne.df.tsv.gz",
    "batch1_genes": f"{BASE}/GSE96583_batch1.genes.tsv.gz",
}


def download(url: str, dest: Path):
    if dest.exists():
        print(f"  already exists: {dest.name}")
        return
    print(f"  downloading {dest.name}...")
    urllib.request.urlretrieve(url, dest)
    print(f"  done ({dest.stat().st_size // 1024} KB)")


for name, url in FILES.items():
    fname = url.split("/")[-1]
    download(url, OUT_DIR / fname)

# Decompress and inspect tsne/metadata files
for batch in [1, 2]:
    gz = OUT_DIR / f"GSE96583_batch{batch}.total.tsne.df.tsv.gz"
    if not gz.exists():
        continue
    tsv = OUT_DIR / f"GSE96583_batch{batch}.total.tsne.df.tsv"
    if not tsv.exists():
        with gzip.open(gz, 'rb') as f_in, open(tsv, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    df = pd.read_csv(tsv, sep="\t")
    print(f"\n=== batch{batch} tsne/metadata ===")
    print(f"Shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    for col in df.columns:
        n = df[col].nunique()
        if n <= 30:
            print(f"  {col} ({n}): {df[col].value_counts().to_dict()}")
        else:
            print(f"  {col}: {n} unique (sample: {df[col].dropna().iloc[:3].tolist()})")
