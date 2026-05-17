"""
Pseudobulk differential expression: vWAT Areg cells, TC vs SC
Mouse-level inference (N=3 SC, N=4 TC) using DESeq2 via pydeseq2.

Replaces cell-level Wilcoxon DE (pseudoreplication) with proper pseudobulk.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings
warnings.filterwarnings('ignore')

DATA_PATH = '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/GSE183288_Single_cell_atlas.h5ad'
OUT_DIR = '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/python_pipeline/outputs'

# ── 1. Load and filter to vWAT Areg cells ────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)

mask = (
    (adata.obs['tissue'] == 'vWAT') &
    (adata.obs['cell_state_label'] == 'Areg')
)
areg = adata[mask].copy()
print(f"vWAT Areg cells: {areg.shape[0]}")
print(areg.obs['intervention_group'].value_counts())

# ── 2. Keep only TC and SC mice ───────────────────────────────────────────────
tc_sc_mask = areg.obs['intervention_group'].isin(['TC', 'SC'])
areg = areg[tc_sc_mask].copy()
print(f"\nTC+SC only: {areg.shape[0]} cells")
print(areg.obs.groupby(['intervention_group', 'sample_name']).size())

# ── 3. Build pseudobulk: sum raw counts per mouse ────────────────────────────
print("\nBuilding pseudobulk...")
samples = areg.obs['sample_name'].unique()
sample_meta = []
count_rows = []

for sample in sorted(samples):
    cells = areg[areg.obs['sample_name'] == sample]
    group = cells.obs['intervention_group'].iloc[0]
    n_cells = cells.shape[0]

    # Use raw counts
    if 'counts' in cells.layers:
        X = cells.layers['counts']
    else:
        X = cells.X

    if issparse(X):
        counts = np.array(X.sum(axis=0)).flatten()
    else:
        counts = X.sum(axis=0).flatten()

    count_rows.append(counts.astype(int))
    sample_meta.append({
        'sample': sample,
        'group': group,
        'n_cells': n_cells
    })

counts_df = pd.DataFrame(
    np.array(count_rows),
    index=[m['sample'] for m in sample_meta],
    columns=areg.var_names
)
meta_df = pd.DataFrame(sample_meta).set_index('sample')

print("\nPseudobulk sample table:")
print(meta_df)

# ── 4. Filter lowly expressed genes ──────────────────────────────────────────
# Keep genes with ≥10 counts in at least 3 samples
min_counts = 10
min_samples = 3
gene_mask = (counts_df >= min_counts).sum(axis=0) >= min_samples
counts_df = counts_df.loc[:, gene_mask]
print(f"\nGenes after filtering (≥{min_counts} counts in ≥{min_samples} samples): {counts_df.shape[1]}")

# ── 5. Run DESeq2: TC vs SC ───────────────────────────────────────────────────
print("\nRunning DESeq2 (TC vs SC)...")

dds = DeseqDataSet(
    counts=counts_df,
    metadata=meta_df,
    design_factors=['group'],
    ref_level=['group', 'SC'],
    quiet=True
)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=['group', 'TC', 'SC'], quiet=True)
stat_res.summary()
results = stat_res.results_df.copy()

# ── 6. Report results ─────────────────────────────────────────────────────────
results_sorted = results.sort_values('padj')
sig_up = results[(results['padj'] < 0.05) & (results['log2FoldChange'] > 0)]
sig_down = results[(results['padj'] < 0.05) & (results['log2FoldChange'] < 0)]

print(f"\n{'='*60}")
print(f"DESeq2 results: TC vs SC, vWAT Areg pseudobulk")
print(f"{'='*60}")
print(f"Total genes tested: {len(results)}")
print(f"Significant upregulated (padj<0.05, logFC>0): {len(sig_up)}")
print(f"Significant downregulated (padj<0.05, logFC<0): {len(sig_down)}")

print(f"\nTop 20 upregulated genes:")
print(sig_up.sort_values('padj').head(20)[['log2FoldChange', 'padj']].to_string())

print(f"\nTop 20 downregulated genes:")
print(sig_down.sort_values('padj').head(20)[['log2FoldChange', 'padj']].to_string())

# ── 7. Check our key genes specifically ──────────────────────────────────────
genes_of_interest = ['Fndc5', 'Ppargc1b', 'Dbp', 'Nr1d1', 'Nr1d2', 'Tef',
                     'Per3', 'Bhlhe41', 'Cdkn1a', 'Clock', 'Arntl']

print(f"\n{'='*60}")
print("Key genes of interest:")
print(f"{'='*60}")
for gene in genes_of_interest:
    if gene in results.index:
        row = results.loc[gene]
        sig = '***' if row['padj'] < 0.001 else ('**' if row['padj'] < 0.01 else ('*' if row['padj'] < 0.05 else 'ns'))
        print(f"  {gene:12s}  logFC={row['log2FoldChange']:+.2f}  padj={row['padj']:.3g}  {sig}")
    else:
        print(f"  {gene:12s}  not expressed / filtered out")

# ── 8. CLOCK target enrichment on pseudobulk results ─────────────────────────
from scipy.stats import hypergeom, fisher_exact

clock_targets = [
    'Dbp', 'Nr1d1', 'Nr1d2', 'Tef', 'Per1', 'Per2', 'Per3',
    'Cry1', 'Cry2', 'Rora', 'Rorc', 'Bhlhe40', 'Bhlhe41',
    'Cdkn1a', 'Nampt', 'Wee1'
]

n_tested = len(results)
n_sig_up = len(sig_up)
n_clock = len(clock_targets)
clock_in_sig = [g for g in clock_targets if g in sig_up.index]
n_clock_in_sig = len(clock_in_sig)

print(f"\n{'='*60}")
print("CLOCK target enrichment (pseudobulk):")
print(f"{'='*60}")
print(f"Background genes tested: {n_tested}")
print(f"Significantly upregulated: {n_sig_up}")
print(f"CLOCK targets in curated set: {n_clock}")
print(f"CLOCK targets in sig upregulated: {n_clock_in_sig}")
print(f"CLOCK targets found: {clock_in_sig}")

if n_sig_up > 0 and n_clock_in_sig > 0:
    expected = n_sig_up * n_clock / n_tested
    fold_enrichment = (n_clock_in_sig / n_sig_up) / (n_clock / n_tested)

    # Hypergeometric p-value
    p_hyper = hypergeom.sf(n_clock_in_sig - 1, n_tested, n_clock, n_sig_up)

    # Fisher's exact
    contingency = [
        [n_clock_in_sig, n_clock - n_clock_in_sig],
        [n_sig_up - n_clock_in_sig, n_tested - n_sig_up - (n_clock - n_clock_in_sig)]
    ]
    _, p_fisher = fisher_exact(contingency, alternative='greater')

    print(f"\nExpected by chance: {expected:.2f}")
    print(f"Fold enrichment: {fold_enrichment:.1f}x")
    print(f"Hypergeometric p: {p_hyper:.3g}")
    print(f"Fisher's exact p: {p_fisher:.3g}")
else:
    print("\nInsufficient significant genes for enrichment test.")
    print("NOTE: This is the pseudoreplication problem in action.")
    print("With N=3 vs N=4 mice, pseudobulk DESeq2 has very low power.")

# ── 9. Save results ───────────────────────────────────────────────────────────
out_path = f"{OUT_DIR}/pseudobulk_de_vwat_areg_tc_vs_sc.csv"
results.to_csv(out_path)
print(f"\nFull results saved to: {out_path}")

print(f"\n{'='*60}")
print("COMPARISON TO CELL-LEVEL WILCOXON:")
print(f"{'='*60}")
print(f"  Wilcoxon significant upregulated: 619 genes")
print(f"  Pseudobulk significant upregulated: {len(sig_up)} genes")
print(f"  Difference reflects pseudoreplication inflation in Wilcoxon")
