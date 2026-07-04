"""
Generate human validation figure from GSE295708 (Miranda et al. 2025).
Shows PPARGC1B and NR1D2 suppressed by obesity and restored by weight loss
in human subcutaneous fat APC (Areg-equivalent) cells.
FNDC5 shown as below-detection control.
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.io
import gzip
import gc
from pathlib import Path
from scipy.sparse import issparse, csr_matrix
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

DATA_DIR = Path('/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/GSE295708_RAW')
OUT_DIR  = Path('outputs')

POOL_META = {
    'Pool_1':  {'gsm': 'GSM8955403', 'group': 'obese'},
    'Pool_2':  {'gsm': 'GSM8955404', 'group': 'lean'},
    'Pool_3':  {'gsm': 'GSM8955405', 'group': 'weight_loss'},
    'Pool_4':  {'gsm': 'GSM8955406', 'group': 'obese'},
    'Pool_5':  {'gsm': 'GSM8955407', 'group': 'lean'},
    'Pool_6':  {'gsm': 'GSM8955408', 'group': 'weight_loss'},
    'Pool_7':  {'gsm': 'GSM8955409', 'group': 'obese'},
    'Pool_8':  {'gsm': 'GSM8955410', 'group': 'lean'},
    'Pool_9':  {'gsm': 'GSM8955411', 'group': 'weight_loss'},
    'Pool_10': {'gsm': 'GSM8955412', 'group': 'obese'},
    'Pool_11': {'gsm': 'GSM8955413', 'group': 'lean'},
    'Pool_12': {'gsm': 'GSM8955414', 'group': 'weight_loss'},
    'Pool_13': {'gsm': 'GSM8955415', 'group': 'obese'},
    'Pool_14': {'gsm': 'GSM8955416', 'group': 'lean'},
    'Pool_15': {'gsm': 'GSM8955417', 'group': 'weight_loss'},
    'Pool_16': {'gsm': 'GSM8955418', 'group': 'obese'},
    'Pool_17': {'gsm': 'GSM8955419', 'group': 'lean'},
    'Pool_18': {'gsm': 'GSM8955420', 'group': 'weight_loss'},
}

# female-only pools for covariate
FEMALE_ONLY_POOLS = {'Pool_3', 'Pool_6', 'Pool_9'}

print('Loading pools...')
adatas = []
for pool_id, meta in POOL_META.items():
    gsm = meta['gsm']
    mtx_file      = DATA_DIR / f'{gsm}_{pool_id}_matrix.mtx.gz'
    barcodes_file = DATA_DIR / f'{gsm}_{pool_id}_barcodes.tsv.gz'
    features_file = DATA_DIR / f'{gsm}_{pool_id}_features.tsv.gz'

    with gzip.open(mtx_file, 'rb') as f:
        mat = csr_matrix(scipy.io.mmread(f), dtype='float32').T
    counts_per_cell = np.asarray(mat.sum(axis=1)).ravel()
    mat = mat[counts_per_cell >= 500]
    gc.collect()

    features = pd.read_csv(features_file, sep='\t', header=None, compression='gzip')
    barcodes = pd.read_csv(barcodes_file, header=None, compression='gzip')[0].values
    barcodes = barcodes[counts_per_cell >= 500]

    adata = sc.AnnData(X=mat)
    adata.var_names = features[1].values
    adata.var['ensembl_id'] = features[0].values
    adata.var_names_make_unique()
    adata.obs_names = [f'{pool_id}_{bc}' for bc in barcodes]
    adata.obs['pool_id'] = pool_id
    adata.obs['group']   = meta['group']
    adatas.append(adata)
    print(f'  {pool_id} ({meta["group"]}): {adata.shape[0]:,} cells')
    del mat; gc.collect()

combined = sc.concat(adatas, join='inner')
del adatas; gc.collect()
print(f'Combined: {combined.shape[0]:,} cells x {combined.shape[1]:,} genes')

# QC
combined.var['mt'] = combined.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], inplace=True)
mask = (
    (combined.obs['n_genes_by_counts'] >= 200) &
    (combined.obs['n_genes_by_counts'] <= 7000) &
    (combined.obs['pct_counts_mt'] <= 5)
)
combined = combined[mask].copy(); gc.collect()
print(f'After QC: {combined.shape[0]:,} cells')

# save raw counts before normalizing
combined.layers['counts'] = combined.X.copy()

# Normalize
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)

# Cluster
sc.pp.highly_variable_genes(combined, n_top_genes=2000, subset=False)
sc.tl.pca(combined, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=20)
sc.tl.leiden(combined, resolution=0.5, key_added='leiden')
print(f'Leiden clusters: {combined.obs["leiden"].nunique()}')

# Find APC cluster (F3-max, PDGFRA+, PTPRC-, PECAM1-)
markers = ['F3', 'PDGFRA', 'PTPRC', 'PECAM1', 'FNDC5', 'PPARGC1B', 'NR1D2', 'NR1D1']
markers = [g for g in markers if g in combined.var_names]

rows = []
for cluster in sorted(combined.obs['leiden'].unique(), key=int):
    sub = combined[combined.obs['leiden'] == cluster]
    row = {'cluster': cluster, 'n_cells': sub.shape[0]}
    for gene in markers:
        vals = sub[:, gene].X
        if issparse(vals): vals = vals.toarray().ravel()
        row[gene] = float(vals.mean())
    rows.append(row)
marker_df = pd.DataFrame(rows).set_index('cluster')

apc_cluster = marker_df['F3'].idxmax()
print(f'\nAPC cluster: {apc_cluster}')
print(marker_df.loc[apc_cluster, markers].round(4))

apc = combined[combined.obs['leiden'] == apc_cluster].copy()
print(f'APC cells: {apc.shape[0]:,}')

# Per-pool means (log-normalized)
key_genes = [g for g in ['PPARGC1B', 'NR1D2', 'NR1D1', 'FNDC5'] if g in apc.var_names]
rows = []
for pool_id in sorted(apc.obs['pool_id'].unique()):
    sub = apc[apc.obs['pool_id'] == pool_id]
    group = POOL_META[pool_id]['group']
    row = {'pool_id': pool_id, 'group': group, 'n_cells': sub.shape[0],
           'female_only': int(pool_id in FEMALE_ONLY_POOLS)}
    for gene in key_genes:
        vals = sub[:, gene].X
        if issparse(vals): vals = vals.toarray().ravel()
        row[gene] = float(vals.mean())
    rows.append(row)
pool_df = pd.DataFrame(rows).set_index('pool_id')

print('\nPer-pool means:')
print(pool_df.round(4).to_string())

# DESeq2
apc_raw = apc.copy()
apc_raw.X = apc_raw.layers['counts']

pseudobulk_rows, pool_ids_ordered = [], []
for pool_id in sorted(apc_raw.obs['pool_id'].unique()):
    sub = apc_raw[apc_raw.obs['pool_id'] == pool_id]
    counts = np.asarray(sub.X.sum(axis=0)).ravel().astype(int)
    pseudobulk_rows.append(counts)
    pool_ids_ordered.append(pool_id)

pb_counts = pd.DataFrame(np.array(pseudobulk_rows),
                          index=pool_ids_ordered, columns=apc_raw.var_names)
pb_meta = pd.DataFrame(
    [{'group': POOL_META[p]['group'],
      'female_only': int(p in FEMALE_ONLY_POOLS),
      'n_cells': pool_df.loc[p, 'n_cells']} for p in pool_ids_ordered],
    index=pool_ids_ordered
)

# drop pools with fewer than 50 APC cells — too sparse for reliable pseudobulk
keep = pb_meta['n_cells'] >= 50
print(f'\nDropping {(~keep).sum()} pools with <50 APC cells:')
print(pb_meta[~keep][['group','n_cells']])
pb_counts = pb_counts[keep]
pb_meta   = pb_meta[keep]

gene_totals = pb_counts.sum(axis=0)
pb_filtered = pb_counts[gene_totals[gene_totals >= 10].index]
print(f'\nGenes tested: {pb_filtered.shape[1]:,}')

results_all = {}
for contrast_name, groups, ref in [
    ('lean_vs_obese',  ['lean', 'obese'],        'obese'),
    ('wl_vs_obese',    ['weight_loss', 'obese'],  'obese'),
]:
    mask = pb_meta['group'].isin(groups)
    counts_sub = pb_filtered.loc[mask]
    meta_sub   = pb_meta.loc[mask].copy()
    meta_sub['group'] = pd.Categorical(meta_sub['group'], categories=groups)

    # only include female_only covariate if it has variance in this subset
    use_female_only = meta_sub['female_only'].nunique() > 1
    design = ['group', 'female_only'] if use_female_only else ['group']
    print(f'\n{contrast_name}: design={design}')

    dds = DeseqDataSet(counts=counts_sub, metadata=meta_sub,
                       design_factors=design,
                       ref_level=['group', ref], quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['group', groups[0], ref],
                       inference=DefaultInference(), quiet=True)
    stats.summary()
    results_all[contrast_name] = stats.results_df.copy()
    print(f'{contrast_name}: {(stats.results_df["padj"]<0.05).sum()} sig genes')

print('\n=== Key genes ===')
for gene in ['PPARGC1B', 'NR1D2', 'NR1D1', 'FNDC5']:
    for contrast_name, res in results_all.items():
        if gene in res.index:
            row = res.loc[gene]
            padj_str = f'{row["padj"]:.2e}' if pd.notna(row['padj']) else 'NaN (below detection)'
            print(f'{gene:12s} {contrast_name:20s} logFC={row["log2FoldChange"]:+.3f}  padj={padj_str}')

# --- Figure ---
colors = {'lean': '#2166AC', 'obese': '#D73027', 'weight_loss': '#969696'}
group_order  = ['lean', 'obese', 'weight_loss']
group_labels = ['Lean', 'Obese', 'Wt. Loss']
group_x      = {g: i for i, g in enumerate(group_order)}

plot_genes = [g for g in ['PPARGC1B', 'NR1D2', 'FNDC5'] if g in key_genes]
fig, axes = plt.subplots(1, len(plot_genes), figsize=(3.2 * len(plot_genes), 4.5))
if len(plot_genes) == 1: axes = [axes]

np.random.seed(42)
for ax, gene in zip(axes, plot_genes):
    for grp in group_order:
        vals = pool_df[pool_df['group'] == grp][gene].values
        jitter = np.random.uniform(-0.1, 0.1, len(vals))
        ax.scatter(np.full(len(vals), group_x[grp]) + jitter, vals,
                   color=colors[grp], s=60, zorder=3,
                   edgecolors='black', linewidths=0.5)
        ax.hlines(vals.mean(), group_x[grp]-0.2, group_x[grp]+0.2,
                  colors=colors[grp], lw=2.5, zorder=4)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(group_labels, fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)

    # annotate padj for lean vs obese and wl vs obese
    y_max = pool_df[gene].max()
    y_range = pool_df[gene].max() - pool_df[gene].min()
    offset = y_range * 0.12

    for contrast_name, x_pos, label in [
        ('lean_vs_obese', 0.5, 'lean vs obese'),
        ('wl_vs_obese',   1.5, 'wl vs obese'),
    ]:
        if gene in results_all.get(contrast_name, pd.DataFrame()).index:
            padj = results_all[contrast_name].loc[gene, 'padj']
            if pd.notna(padj) and padj < 0.05:
                sig_str = f'padj={padj:.2e}'
                ax.annotate(sig_str, xy=(x_pos, y_max + offset),
                            ha='center', fontsize=7, color='#333333')
            elif not pd.notna(padj):
                ax.annotate('n.d.', xy=(x_pos, y_max + offset),
                            ha='center', fontsize=7, color='#999999',
                            style='italic')

    if gene == 'FNDC5':
        ax.set_title('FNDC5\n(not regulated, n.s.)', fontsize=9, color='#888888')
    else:
        ax.set_title(gene, fontsize=10, fontweight='bold')

    if ax == axes[0]:
        ax.set_ylabel('Mean log-normalized expression\n(APC cluster)', fontsize=8)

fig.suptitle('Human subcutaneous fat progenitors (GSE295708, Miranda et al. 2025)\n'
             '70 donors · 18 pools · lean / obese / post-bariatric weight loss',
             fontsize=9, y=1.02)
plt.tight_layout()
out_path = OUT_DIR / 'fig_human_validation.pdf'
plt.savefig(out_path, bbox_inches='tight')
print(f'\nSaved to {out_path}')
