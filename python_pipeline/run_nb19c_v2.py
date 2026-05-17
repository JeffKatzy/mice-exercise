"""
NB19c v2 — Proper Areg-equivalent cell selection + pseudobulk DESeq2

Fixes vs v1:
1. Cell selection: F3>0 AND PDGFRA>0 AND PTPRC=0 AND PECAM1=0 (Areg-equivalent, not whole cluster)
2. Statistics: pseudobulk DESeq2 (pydeseq2) with sex covariate, not group means
3. Sex covariate: Pool_14 (Lean Rep5), Pool_16 (Obese Rep6), Pool_18 (WL Rep6) are female-only pools
"""
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

# Sex metadata from series matrix file (female-only pools confirmed from GEO)
POOL_META = {
    'Pool_1':  {'gsm': 'GSM8955403', 'group': 'obese',       'female_only': False},
    'Pool_2':  {'gsm': 'GSM8955404', 'group': 'lean',        'female_only': False},
    'Pool_3':  {'gsm': 'GSM8955405', 'group': 'weight_loss', 'female_only': False},
    'Pool_4':  {'gsm': 'GSM8955406', 'group': 'obese',       'female_only': False},
    'Pool_5':  {'gsm': 'GSM8955407', 'group': 'lean',        'female_only': False},
    'Pool_6':  {'gsm': 'GSM8955408', 'group': 'weight_loss', 'female_only': False},
    'Pool_7':  {'gsm': 'GSM8955409', 'group': 'obese',       'female_only': False},
    'Pool_8':  {'gsm': 'GSM8955410', 'group': 'lean',        'female_only': False},
    'Pool_9':  {'gsm': 'GSM8955411', 'group': 'weight_loss', 'female_only': False},
    'Pool_10': {'gsm': 'GSM8955412', 'group': 'obese',       'female_only': False},
    'Pool_11': {'gsm': 'GSM8955413', 'group': 'lean',        'female_only': False},
    'Pool_12': {'gsm': 'GSM8955414', 'group': 'weight_loss', 'female_only': False},
    'Pool_13': {'gsm': 'GSM8955415', 'group': 'obese',       'female_only': False},
    'Pool_14': {'gsm': 'GSM8955416', 'group': 'lean',        'female_only': True},   # Female only
    'Pool_15': {'gsm': 'GSM8955417', 'group': 'weight_loss', 'female_only': False},
    'Pool_16': {'gsm': 'GSM8955418', 'group': 'obese',       'female_only': True},   # Female only
    'Pool_17': {'gsm': 'GSM8955419', 'group': 'lean',        'female_only': False},
    'Pool_18': {'gsm': 'GSM8955420', 'group': 'weight_loss', 'female_only': True},   # Female only
}

# ── STEP 1: LOAD ──────────────────────────────────────────────────────────────
print('Loading pools...')
adatas = []
for pool_id, meta in POOL_META.items():
    gsm = meta['gsm']
    mtx_f = DATA_DIR / f'{gsm}_{pool_id}_matrix.mtx.gz'
    bar_f = DATA_DIR / f'{gsm}_{pool_id}_barcodes.tsv.gz'
    fea_f = DATA_DIR / f'{gsm}_{pool_id}_features.tsv.gz'

    with gzip.open(mtx_f, 'rb') as f:
        mat = csr_matrix(scipy.io.mmread(f), dtype='float32').T

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
    ad.obs['pool_id']     = pool_id
    ad.obs['group']       = meta['group']
    ad.obs['female_only'] = meta['female_only']
    adatas.append(ad)
    print(f'  {pool_id} ({meta["group"]:12s}): {ad.shape[0]:,} cells')
    del mat; gc.collect()

combined = sc.concat(adatas, join='inner')
del adatas; gc.collect()
print(f'Combined: {combined.shape[0]:,} cells x {combined.shape[1]:,} genes')

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

# ── STEP 3: NORMALIZE + CLUSTER ──────────────────────────────────────────────
print('\nNormalizing and clustering...')
combined.layers['counts'] = combined.X.copy()  # save raw counts for pseudobulk
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
combined.layers['lognorm'] = combined.X.copy()

sc.pp.highly_variable_genes(combined, n_top_genes=2000, subset=False)
gc.collect()
sc.tl.pca(combined, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=20)
sc.tl.leiden(combined, resolution=0.5, key_added='leiden')
print(f'Clusters: {combined.obs["leiden"].nunique()}')

# ── STEP 4: IDENTIFY PROGENITOR CLUSTER ──────────────────────────────────────
print('\nIdentifying progenitor cluster...')
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
show = [c for c in ['F3','PDGFRA','FNDC5','ITGB5','PTPRC','PECAM1','n_cells'] if c in marker_df.columns]
print(marker_df[show].sort_values('F3', ascending=False).round(3).to_string())

progenitor_cluster = marker_df['F3'].idxmax()
prog = combined[combined.obs['leiden'] == progenitor_cluster].copy()
print(f'\nProgenitor cluster: {progenitor_cluster} ({prog.shape[0]:,} cells)')

# ── STEP 5: SELECT AREG-EQUIVALENT CELLS ─────────────────────────────────────
print('\nSelecting Areg-equivalent cells (F3>0, PDGFRA>0, PTPRC=0, PECAM1=0)...')

def get_expr(adata, gene):
    if gene not in adata.var_names:
        return np.zeros(adata.shape[0])
    v = adata[:, gene].X
    if issparse(v): v = v.toarray().ravel()
    return v

f3     = get_expr(prog, 'F3')
pdgfra = get_expr(prog, 'PDGFRA')
ptprc  = get_expr(prog, 'PTPRC')
pecam1 = get_expr(prog, 'PECAM1')

# Areg-equivalent: F3-positive, PDGFRA-positive, no immune, no endothelial
areg_mask = (f3 > 0) & (pdgfra > 0) & (ptprc == 0) & (pecam1 == 0)
areg = prog[areg_mask].copy()

print(f'Progenitor cluster total: {prog.shape[0]:,}')
print(f'Areg-equivalent (filtered): {areg.shape[0]:,}  ({100*areg_mask.mean():.1f}% of cluster)')
print('\nAreg cells per pool:')
pool_counts = areg.obs.groupby(['pool_id','group']).size().sort_index()
print(pool_counts.to_string())

# ── STEP 6: GROUP MEANS (for reference) ──────────────────────────────────────
print('\nGroup means in Areg-equivalent cells:')
key_genes = [g for g in ['FNDC5','ITGB5','PPARGC1B','NR1D1','F3'] if g in areg.var_names]

rows = []
for pool_id in sorted(areg.obs['pool_id'].unique()):
    sub = areg[areg.obs['pool_id']==pool_id]
    meta = POOL_META[pool_id]
    row = {'pool_id': pool_id, 'group': meta['group'],
           'female_only': meta['female_only'], 'n_cells': sub.shape[0]}
    for g in key_genes:
        v = sub[:,g].X
        if issparse(v): v = v.toarray().ravel()
        row[g] = float(v.mean())
    rows.append(row)

pool_df = pd.DataFrame(rows).set_index('pool_id')
print(pool_df.round(4).to_string())
print('\nGroup means:')
print(pool_df.groupby('group')[key_genes].mean().round(4).to_string())

# ── STEP 7: PSEUDOBULK DESEQ2 ────────────────────────────────────────────────
print('\n── Pseudobulk DESeq2 ──')
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    HAVE_DESEQ2 = True
except ImportError:
    HAVE_DESEQ2 = False
    print('pydeseq2 not installed. Run: pip install pydeseq2')

if HAVE_DESEQ2:
    # Sum raw counts per pool in Areg cells
    areg.X = areg.layers['counts']  # use raw counts

    pb_rows, pb_meta_rows, pool_ids_ord = [], [], []
    for pool_id in sorted(areg.obs['pool_id'].unique()):
        sub = areg[areg.obs['pool_id']==pool_id]
        if sub.shape[0] < 10:
            print(f'  Skipping {pool_id} — only {sub.shape[0]} Areg cells')
            continue
        counts = np.asarray(sub.X.sum(axis=0)).ravel().astype(int)
        pb_rows.append(counts)
        meta = POOL_META[pool_id]
        pb_meta_rows.append({
            'group': meta['group'],
            'female_only': 'yes' if meta['female_only'] else 'no'
        })
        pool_ids_ord.append(pool_id)
        print(f'  {pool_id} ({meta["group"]:12s}): {sub.shape[0]:,} Areg cells, {int(counts.sum()):,} total counts')

    pb_counts = pd.DataFrame(np.array(pb_rows), index=pool_ids_ord, columns=areg.var_names)
    pb_meta   = pd.DataFrame(pb_meta_rows, index=pool_ids_ord)

    # Filter low-count genes
    keep = pb_counts.sum(axis=0) >= 10
    pb_counts = pb_counts.loc[:, keep]
    print(f'\nGenes kept (total>=10 across pools): {keep.sum():,}')
    print(f'Pools in pseudobulk: {len(pool_ids_ord)}')
    print(f'Metadata:\n{pb_meta.to_string()}')

    results = {}
    for contrast_name, groups, ref in [
        ('lean_vs_obese',   ['lean','obese'],       'obese'),
        ('wl_vs_obese',     ['weight_loss','obese'], 'obese'),
        ('lean_vs_wl',      ['lean','weight_loss'],  'weight_loss'),
    ]:
        mask = pb_meta['group'].isin(groups)
        if mask.sum() < 4:
            print(f'\nSkipping {contrast_name} — not enough pools')
            continue
        counts_sub = pb_counts.loc[mask]
        meta_sub   = pb_meta.loc[mask].copy()

        try:
            dds = DeseqDataSet(
                counts=counts_sub,
                metadata=meta_sub,
                design_factors=['group', 'female_only'],
                ref_level=['group', ref],
                quiet=True,
            )
            dds.deseq2()
            stats = DeseqStats(
                dds,
                contrast=['group', groups[0], ref],
                inference=DefaultInference(),
                quiet=True,
            )
            stats.summary()
            results[contrast_name] = stats.results_df.copy()
            n_sig = (stats.results_df['padj'] < 0.05).sum()
            print(f'\n{contrast_name}: {n_sig} significant genes (padj<0.05)')
        except Exception as e:
            print(f'\n{contrast_name} failed: {e}')

    # Report key genes
    if results:
        genes_of_interest = ['FNDC5','ITGB5','PPARGC1B','NR1D1','NR1D2','F3','PDGFRA']
        print('\n=== Key genes across contrasts ===')
        print(f'{"Gene":<12} {"Contrast":<22} {"logFC":>7}  {"padj":>10}  {"baseMean":>10}')
        print('-' * 65)
        for gene in genes_of_interest:
            for cname, res in results.items():
                if gene in res.index:
                    row = res.loc[gene]
                    lfc  = row.get('log2FoldChange', float('nan'))
                    padj = row.get('padj', float('nan'))
                    bm   = row.get('baseMean', float('nan'))
                    sig  = ' *' if (not np.isnan(padj) and padj < 0.05) else ''
                    print(f'{gene:<12} {cname:<22} {lfc:+7.3f}  {padj:10.3g}  {bm:10.1f}{sig}')

# ── STEP 8: STRIP PLOT ────────────────────────────────────────────────────────
colors = {'lean':'#2166AC','obese':'#D73027','weight_loss':'#F4A582'}
grp_x  = {'lean':0,'obese':1,'weight_loss':2}

fig, axes = plt.subplots(1, len(key_genes), figsize=(3.5*len(key_genes), 5))
if len(key_genes)==1: axes=[axes]

for ax, gene in zip(axes, key_genes):
    for grp in ['lean','obese','weight_loss']:
        sub = pool_df[pool_df['group']==grp]
        vals = sub[gene].values
        x = [grp_x[grp]]*len(vals)
        ax.scatter(x, vals, color=colors[grp], s=70, zorder=3)
        ax.hlines(vals.mean(), grp_x[grp]-0.18, grp_x[grp]+0.18, colors=colors[grp], lw=2.5)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(['Lean','Obese','Wt Loss'], fontsize=8, rotation=15)
    ax.set_title(gene, fontsize=10)
    ax.spines[['top','right']].set_visible(False)

axes[0].set_ylabel('Mean log-norm expression (Areg-equivalent cells)')
plt.suptitle('GSE295708 — Human SAT Areg-equivalent cells\n'
             'F3>0, PDGFRA>0, PTPRC=0, PECAM1=0. Each dot = one pool.', fontsize=9)
plt.tight_layout()
plt.savefig(OUT_DIR / 'areg_three_groups_v2.png', dpi=150, bbox_inches='tight')
print('\nSaved areg_three_groups_v2.png')
print('\nDone.')
