"""
Download Mathys et al. 2019 metadata only and print schema.
Run once before building the full download pipeline.
"""

import synapseclient
import pandas as pd
from pathlib import Path

OUT_DIR = Path(__file__).parent / "dataset" / "mathys"
OUT_DIR.mkdir(parents=True, exist_ok=True)

syn = synapseclient.login()

print("Downloading filtered_column_metadata.txt (syn18686372)...")
f = syn.get("syn18686372", downloadLocation=str(OUT_DIR))
meta_path = Path(f.path)
print(f"Saved to: {meta_path}\n")

meta = pd.read_csv(meta_path, sep="\t")
print(f"Shape: {meta.shape}")
print(f"\nColumns:\n{meta.columns.tolist()}")
print(f"\nDtypes:\n{meta.dtypes.to_string()}")

for col in meta.columns:
    n_unique = meta[col].nunique()
    if n_unique <= 20:
        print(f"\n{col} ({n_unique} unique):\n  {meta[col].value_counts().to_dict()}")
    else:
        print(f"\n{col}: {n_unique} unique values (sample: {meta[col].dropna().iloc[:5].tolist()})")
