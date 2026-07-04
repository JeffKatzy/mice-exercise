import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.sparse import issparse
from pathlib import Path

DATA_PATH = Path('../GSE183288_Single_cell_atlas.h5ad')
OUT_DIR   = Path('outputs')

CONDITION_COLORS = {'SC': '#2166AC', 'TC': '#4DAF4A', 'SH': '#E41A1C', 'TH': '#FF7F00'}
CONDITION_ORDER  = ['SC', 'TC']
GENE             = 'Fndc5'

adata = sc.read_h5ad(DATA_PATH)

def per_mouse_means(adata, tissue, cell_state, conditions, gene):
    mask = (
        (adata.obs['tissue'] == tissue) &
        (adata.obs['cell_state_label'] == cell_state) &
        (adata.obs['intervention_group'].isin(conditions))
    )
    sub = adata[mask]
    gene_idx = list(sub.var_names).index(gene)
    X = sub.X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame({
        'expr':      X[:, gene_idx],
        'sample_ID': sub.obs['sample_ID'].values,
        'condition': sub.obs['intervention_group'].values,
    })
    means = df.groupby(['sample_ID', 'condition'])['expr'].mean().reset_index()
    means.columns = ['sample_ID', 'condition', 'mean_expr']
    return means

vwat_means  = per_mouse_means(adata, 'vWAT',  'Areg', CONDITION_ORDER, GENE)
scwat_means = per_mouse_means(adata, 'scWAT', 'Areg', CONDITION_ORDER, GENE)

# --- figure ---
fig, axes = plt.subplots(1, 2, figsize=(6, 4), sharey=False)

np.random.seed(42)

for ax, means, title, show_pval in [
    (axes[0], vwat_means,  'Visceral fat\n(vWAT)',       True),
    (axes[1], scwat_means, 'Subcutaneous fat\n(scWAT)', False),
]:
    for i, cond in enumerate(CONDITION_ORDER):
        vals   = means[means['condition'] == cond]['mean_expr'].values
        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter, vals,
            color=CONDITION_COLORS[cond], s=70,
            edgecolors='black', linewidths=0.6, zorder=3
        )
        ax.hlines(
            np.mean(vals), i - 0.18, i + 0.18,
            colors=CONDITION_COLORS[cond], linewidths=2.5, zorder=4
        )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Sedentary\n(SC)', 'Trained\n(TC)'], fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_ylabel('Fndc5 mean expression\n(log-normalized)' if ax == axes[0] else '', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if show_pval:
        ax.annotate('p = 0.029\n(permutation)', xy=(0.5, 0.92),
                    xycoords='axes fraction', ha='center', fontsize=8,
                    color='#333333')
    else:
        # show cell count caveat for scWAT
        sc_n = (adata[(adata.obs['tissue']=='scWAT') &
                      (adata.obs['cell_state_label']=='Areg') &
                      (adata.obs['intervention_group'].isin(['SC','TC']))]).n_obs
        ax.annotate(f'n.s.\n({sc_n} cells total)', xy=(0.5, 0.92),
                    xycoords='axes fraction', ha='center', fontsize=8,
                    color='#888888')

fig.suptitle('Fndc5 in Areg cells: visceral vs subcutaneous fat',
             fontsize=11, fontweight='bold', y=1.02)

plt.tight_layout()
out_path = OUT_DIR / 'fig_depot_specificity.pdf'
plt.savefig(out_path, bbox_inches='tight')
print(f'Saved to {out_path}')

# print the per-mouse numbers for reference
print('\nvWAT Areg per mouse:')
print(vwat_means.sort_values(['condition','mean_expr']).to_string(index=False))
print('\nscWAT Areg per mouse:')
print(scwat_means.sort_values(['condition','mean_expr']).to_string(index=False))
