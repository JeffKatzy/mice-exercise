import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.io
import gzip
import gc
from pathlib import Path
from scipy.sparse import issparse, csr_matrix
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path('/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/GSE295708_RAW')
OUT_DIR  = Path('outputs/nb19c_figures')
OUT_DIR.mkdir(parents=True, exist_ok=True)

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

# ── STEP 1: LOAD ─────────────────────────────────────────────────────────────
print('Loading pools...')
adatas = []
for pool_id, meta in POOL_META.items():
    gsm = meta['gsm']
    mtx_f = DATA_DIR / f'{gsm}_{pool_id}_matrix.mtx.gz'
    bar_f = DATA_DIR / f'{gsm}_{pool_id}_barcodes.tsv.gz'
    fea_f = DATA_DIR / f'{gsm}_{pool_id}_features.tsv.gz'

    with gzip.open(mtx_f, 'rb') as f:
        mat = scipy.io.mmread(f)
        mat = csr_matrix(mat, dtype='float32').T  # (cells x genes), float32 immediately

    counts_per_cell = np.asarray(mat.sum(axis=1)).ravel()
    mask = counts_per_cell >= 500
    mat = mat[mask]
    gc.collect()

    features = pd.read_csv(fea_f, sep='\t', header=None, compression='gzip')
    barcodes = pd.read_csv(bar_f, header=None, compression='gzip')[0].values[mask]

    ad = sc.AnnData(X=mat)
    ad.var_names = features[1].values
    ad.var_names_make_unique()
    ad.obs_names = [f'{pool_id}_{bc}' for bc in barcodes]
    ad.obs['pool_id'] = pool_id
    ad.obs['group']   = meta['group']
    adatas.append(ad)
    print(f'  {pool_id} ({meta["group"]:12s}): {ad.shape[0]:,} cells', flush=True)
    del mat; gc.collect()

print('Concatenating...')
combined = sc.concat(adatas, join='inner')
del adatas; gc.collect()
print(f'Combined: {combined.shape[0]:,} cells x {combined.shape[1]:,} genes')
print(combined.obs.groupby('group').size().to_string())

# ── STEP 2: QC ────────────────────────────────────────────────────────────────
print('\nQC filtering...')
combined.var['mt'] = combined.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], inplace=True)
mask = (
    (combined.obs['n_genes_by_counts'] >= 200) &
    (combined.obs['n_genes_by_counts'] <= 7000) &
    (combined.obs['pct_counts_mt'] <= 5)
)
combined = combined[mask].copy()
gc.collect()
print(f'After QC: {combined.shape[0]:,} cells')
print(combined.obs.groupby('group').size().to_string())

# ── STEP 3: NORMALIZE + CLUSTER ──────────────────────────────────────────────
print('\nNormalizing...')
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
# Skip scale — would densify. PCA on log-norm sparse is sufficient for clustering.
sc.pp.highly_variable_genes(combined, n_top_genes=2000, subset=False)
gc.collect()
sc.tl.pca(combined, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=20)
print('Clustering...')
sc.tl.leiden(combined, resolution=0.5, key_added='leiden')
print(f'Clusters: {combined.obs["leiden"].nunique()}')

# ── STEP 4: IDENTIFY APC CLUSTER ─────────────────────────────────────────────
print('\nIdentifying APC cluster...')
markers = ['F3','PDGFRA','CD34','FNDC5','ITGB5','PPARGC1B','NR1D1','PTPRC','PECAM1']
markers = [g for g in markers if g in combined.var_names]

rows = []
for cl in sorted(combined.obs['leiden'].unique(), key=int):
    sub = combined[combined.obs['leiden']==cl]
    row = {'cluster': cl, 'n_cells': sub.shape[0]}
    for g in markers:
        v = sub[:,g].X
        if issparse(v): v = v.toarray().ravel()
        row[g] = float(v.mean())
    rows.append(row)

marker_df = pd.DataFrame(rows).set_index('cluster')
show_cols = [c for c in ['F3','PDGFRA','FNDC5','ITGB5','PTPRC','PECAM1','n_cells'] if c in marker_df.columns]
print(marker_df[show_cols].sort_values('F3', ascending=False).round(3).to_string())

apc_cluster = marker_df['F3'].idxmax()
print(f'\nAPC cluster: {apc_cluster}')
apc = combined[combined.obs['leiden']==apc_cluster].copy()
print(f'APC cells: {apc.shape[0]:,}')
print(apc.obs.groupby(['pool_id','group']).size().sort_index().to_string())

# ── STEP 5: PER-POOL FNDC5 ────────────────────────────────────────────────────
print('\nFNDC5 per pool in APC cluster:')
key_genes = [g for g in ['FNDC5','ITGB5','PPARGC1B','NR1D1','F3'] if g in apc.var_names]
rows = []
for pool_id in sorted(apc.obs['pool_id'].unique()):
    sub = apc[apc.obs['pool_id']==pool_id]
    group = POOL_META[pool_id]['group']
    row = {'pool_id': pool_id, 'group': group, 'n_cells': sub.shape[0]}
    for g in key_genes:
        v = sub[:,g].X
        if issparse(v): v = v.toarray().ravel()
        row[g] = float(v.mean())
    rows.append(row)

pool_df = pd.DataFrame(rows).set_index('pool_id')
print(pool_df.round(4).to_string())
print('\nGroup means:')
print(pool_df.groupby('group')[key_genes].mean().round(4).to_string())

lean  = pool_df[pool_df['group']=='lean']['FNDC5'].values
obese = pool_df[pool_df['group']=='obese']['FNDC5'].values
wl    = pool_df[pool_df['group']=='weight_loss']['FNDC5'].values
print(f'\nFNDC5 lean:        {lean.mean():.4f}  (range {lean.min():.4f}-{lean.max():.4f})')
print(f'FNDC5 obese:       {obese.mean():.4f}  (range {obese.min():.4f}-{obese.max():.4f})')
print(f'FNDC5 weight_loss: {wl.mean():.4f}  (range {wl.min():.4f}-{wl.max():.4f})')
print(f'Lean > Obese: {lean.mean() > obese.mean()}')
print(f'WL between lean and obese: {obese.mean() < wl.mean() < lean.mean()}')

# ── STEP 6: PLOT ──────────────────────────────────────────────────────────────
colors = {'lean':'#2166AC','obese':'#D73027','weight_loss':'#F4A582'}
grp_x  = {'lean':0,'obese':1,'weight_loss':2}
fig, axes = plt.subplots(1, len(key_genes), figsize=(3.5*len(key_genes), 5))
if len(key_genes)==1: axes=[axes]
for ax, gene in zip(axes, key_genes):
    for grp in ['lean','obese','weight_loss']:
        vals = pool_df[pool_df['group']==grp][gene].values
        x = [grp_x[grp]]*len(vals)
        ax.scatter(x, vals, color=colors[grp], s=70, zorder=3)
        ax.hlines(vals.mean(), grp_x[grp]-0.18, grp_x[grp]+0.18, colors=colors[grp], lw=2.5)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(['Lean','Obese','Wt Loss'], fontsize=8, rotation=15)
    ax.set_title(gene, fontsize=10)
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('Mean log-normalized expression (APC cluster)')
plt.suptitle('GSE295708 — Human subcutaneous fat APC cells\nEach dot = one pool. Bar = group mean.', fontsize=10)
plt.tight_layout()
plt.savefig(OUT_DIR / 'fndc5_three_groups.png', dpi=150, bbox_inches='tight')
print('Saved fndc5_three_groups.png')
print('\nDone.')
