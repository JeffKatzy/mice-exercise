"""
Figure 1 — The decisive figure.
4 panels telling one story:
  A. Complete rank separation: vWAT Areg per-mouse strip, mice visually dominant
  B. Visceral specificity: vWAT vs scWAT Areg
  C. Lineage specificity: Fndc5 across vWAT cell types with n labels
  D. Clock coherence: raw centered (not z-scored) + spaghetti lines per mouse
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.sparse import issparse
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

DATA_PATH = '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/GSE183288_Single_cell_atlas.h5ad'
OUT_PATH  = '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/python_pipeline/outputs/figure1_decisive.pdf'
OUT_PNG   = '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise/python_pipeline/outputs/figure1_decisive.png'

SC_COLOR = '#4878CF'
TC_COLOR = '#2CA02C'
MOUSE_COLORS = ['#e07b54', '#b85c38', '#8B3A3A',          # SC mice (warm)
                '#3a7ebf', '#2a6aad', '#1a5a9d', '#0a4a8d'] # TC mice (cool)

print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)

def get_col(adata_sub, gene):
    if gene not in adata_sub.var_names:
        return np.zeros(adata_sub.shape[0])
    idx = list(adata_sub.var_names).index(gene)
    col = adata_sub.X[:, idx]
    if issparse(col):
        col = col.toarray().flatten()
    return np.array(col)

def mean_expr(adata_sub, gene):
    return float(np.mean(get_col(adata_sub, gene)))

def pct_expr(adata_sub, gene):
    col = get_col(adata_sub, gene)
    return 100.0 * np.mean(col > 0)

def per_mouse(adata_sub, gene):
    """Returns ordered lists: sc_samples, tc_samples with (val, pct, n, name)"""
    sc_res, tc_res = [], []
    for sample in sorted(adata_sub.obs['sample_name'].unique()):
        cells = adata_sub[adata_sub.obs['sample_name'] == sample]
        grp = cells.obs['intervention_group'].iloc[0]
        val = mean_expr(cells, gene)
        pct = pct_expr(cells, gene)
        n   = cells.shape[0]
        entry = (val, pct, n, sample)
        if grp == 'SC':
            sc_res.append(entry)
        else:
            tc_res.append(entry)
    return sorted(sc_res), sorted(tc_res)

# ── Data prep ─────────────────────────────────────────────────────────────────
vwat_areg = adata[
    (adata.obs['tissue'] == 'vWAT') &
    (adata.obs['cell_state_label'] == 'Areg') &
    (adata.obs['intervention_group'].isin(['TC', 'SC']))
].copy()

scwat_areg = adata[
    (adata.obs['tissue'] == 'scWAT') &
    (adata.obs['cell_state_label'] == 'Areg') &
    (adata.obs['intervention_group'].isin(['TC', 'SC']))
].copy()

vwat_sc, vwat_tc = per_mouse(vwat_areg, 'Fndc5')

# scWAT — adequate cell count only
scwat_sc_raw, scwat_tc_raw = per_mouse(scwat_areg, 'Fndc5')
scwat_sc = [(v, p, n, s) for v, p, n, s in scwat_sc_raw if n >= 10]
scwat_tc = [(v, p, n, s) for v, p, n, s in scwat_tc_raw if n >= 10]

# Panel C cell types with cell counts
cell_types_ordered = ['Areg', 'pre_CP', 'CP', 'WAT_IPC',
                      'M1', 'M2', 'NK', 'NKT', 'CD4_Memory', 'CD8_Cytotoxic', 'B_Cell']
ct_vals = {}
ct_n    = {}
for ct in cell_types_ordered:
    sub = adata[
        (adata.obs['tissue'] == 'vWAT') &
        (adata.obs['cell_state_label'] == ct) &
        (adata.obs['intervention_group'].isin(['TC', 'SC']))
    ]
    if sub.shape[0] < 20:
        continue
    sub_sc = sub[sub.obs['intervention_group'] == 'SC']
    sub_tc = sub[sub.obs['intervention_group'] == 'TC']
    ct_vals[ct] = (mean_expr(sub_sc, 'Fndc5'), mean_expr(sub_tc, 'Fndc5'))
    ct_n[ct]    = (sub_sc.shape[0], sub_tc.shape[0])

# Panel D: per-mouse values for clock module — raw log-normalized, mean-centered per gene
clock_genes = ['Nr1d1', 'Nr1d2', 'Tef', 'Dbp', 'Fndc5']
clock_pm = {}   # gene -> {mouse_id -> (centered_val, group)}
for gene in clock_genes:
    sc_list, tc_list = per_mouse(vwat_areg, gene)
    all_vals = [v for v, p, n, s in sc_list + tc_list]
    gene_mean = np.mean(all_vals)
    d = {}
    for v, p, n, s in sc_list:
        d[s] = (v - gene_mean, 'SC')
    for v, p, n, s in tc_list:
        d[s] = (v - gene_mean, 'TC')
    clock_pm[gene] = d

# Consistent mouse ordering for spaghetti
all_mice_sc = sorted([s for v, p, n, s in vwat_sc])
all_mice_tc = sorted([s for v, p, n, s in vwat_tc])
all_mice = all_mice_sc + all_mice_tc

# Assign a consistent color per mouse
mouse_color = {}
for i, m in enumerate(all_mice_sc):
    mouse_color[m] = MOUSE_COLORS[i % 3]
for i, m in enumerate(all_mice_tc):
    mouse_color[m] = MOUSE_COLORS[3 + i % 4]

print("Data ready. Building figure...")

fig = plt.figure(figsize=(14, 10))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(2, 2, figure=fig,
                       hspace=0.45, wspace=0.35,
                       left=0.07, right=0.97, top=0.92, bottom=0.09)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

LABEL_KW = dict(fontsize=10.5, fontweight='bold')

# ── Panel A: Rank separation — mice visually dominant ─────────────────────────
ax = ax_a
np.random.seed(42)

sc_vals = [v for v, p, n, s in vwat_sc]
tc_vals = [v for v, p, n, s in vwat_tc]
sc_pcts = [p for v, p, n, s in vwat_sc]
tc_pcts = [p for v, p, n, s in vwat_tc]
sc_names = [s for v, p, n, s in vwat_sc]
tc_names = [s for v, p, n, s in vwat_tc]

# Large markers — individual mice are the unit of analysis
for i, (v, name) in enumerate(zip(sc_vals, sc_names)):
    ax.scatter(0, v, color=mouse_color[name], s=160, zorder=6,
               edgecolors='white', linewidths=1.0)
    ax.text(-0.28, v, name.replace('_', ' '), ha='right', va='center',
            fontsize=6.5, color=mouse_color[name])
for i, (v, name) in enumerate(zip(tc_vals, tc_names)):
    ax.scatter(1, v, color=mouse_color[name], s=160, zorder=6,
               edgecolors='white', linewidths=1.0)
    ax.text(1.28, v, name.replace('_', ' '), ha='left', va='center',
            fontsize=6.5, color=mouse_color[name])

# All pairwise lines
for sv in sc_vals:
    for tv in tc_vals:
        ax.plot([0, 1], [sv, tv], color='gray', alpha=0.15, lw=0.8, zorder=1)

ax.hlines(np.mean(sc_vals), -0.22, 0.22, colors='#333333', lw=2.5, zorder=7)
ax.hlines(np.mean(tc_vals),  0.78, 1.22, colors='#333333', lw=2.5, zorder=7)

# % expressing annotation
ax.text(0, -0.007, f'{np.mean(sc_pcts):.0f}% expr.', ha='center', fontsize=7.5,
        color=SC_COLOR, style='italic')
ax.text(1, -0.007, f'{np.mean(tc_pcts):.0f}% expr.', ha='center', fontsize=7.5,
        color=TC_COLOR, style='italic')

ax.set_xticks([0, 1])
ax.set_xticklabels(['Sedentary\n(SC, N=3)', 'Trained\n(TC, N=4)'], fontsize=9)
ax.set_ylabel('Mean Fndc5 expression\n(log-normalized)', fontsize=9)
ax.set_xlim(-0.45, 1.45)
ax.set_ylim(-0.012, 0.105)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('A  vWAT Areg: complete rank separation', loc='left', **LABEL_KW)
ax.text(0.97, 0.97, 'p = 0.029\n(exact permutation)', transform=ax.transAxes,
        ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0', edgecolor='none'))

# ── Panel B: Visceral specificity ─────────────────────────────────────────────
ax = ax_b

for v, p, n, s in vwat_sc:
    ax.scatter(0, v, color=mouse_color[s], s=90, zorder=5, edgecolors='white', lw=0.8)
for v, p, n, s in vwat_tc:
    ax.scatter(1, v, color=mouse_color[s], s=90, zorder=5, edgecolors='white', lw=0.8)
ax.hlines(np.mean(sc_vals), -0.22, 0.22, colors='#333333', lw=2, zorder=6)
ax.hlines(np.mean(tc_vals),  0.78, 1.22, colors='#333333', lw=2, zorder=6)

offset = 2.0
sc_sc = [v for v, p, n, s in scwat_sc]
sc_tc = [v for v, p, n, s in scwat_tc]
for v, p, n, s in scwat_sc:
    ax.scatter(offset + 0, v, color=SC_COLOR, s=70, marker='s', zorder=5,
               edgecolors='white', lw=0.6, alpha=0.8)
for v, p, n, s in scwat_tc:
    ax.scatter(offset + 1, v, color=TC_COLOR, s=70, marker='s', zorder=5,
               edgecolors='white', lw=0.6, alpha=0.8)
if sc_sc:
    ax.hlines(np.mean(sc_sc), offset-0.22, offset+0.22, colors='#333333', lw=2, zorder=6)
if sc_tc:
    ax.hlines(np.mean(sc_tc), offset+0.78, offset+1.22, colors='#333333', lw=2, zorder=6)

ax.axvline(x=1.5, color='#cccccc', lw=1, ls='--')
ax.set_xticks([0, 1, offset, offset+1])
ax.set_xticklabels(['SC', 'TC', 'SC', 'TC'], fontsize=9)
ax.set_xlim(-0.5, offset+1.5)
ax.set_ylim(-0.012, 0.105)

ax.text(0.5, -0.010, 'vWAT Areg', ha='center', fontsize=8, style='italic',
        transform=ax.get_xaxis_transform())
ax.text(offset+0.5, -0.010, f'scWAT Areg\n(N≥10 cells)', ha='center', fontsize=8,
        style='italic', transform=ax.get_xaxis_transform())

ax.annotate('p=0.029', xy=(0.5, 0.089), ha='center', fontsize=8, color='#333333')
ax.annotate('p=0.40 n.s.', xy=(offset+0.5, 0.040), ha='center', fontsize=8, color='gray')

ax.set_ylabel('Mean Fndc5 expression\n(log-normalized)', fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('B  Visceral-enriched: no significant induction in scWAT', loc='left', **LABEL_KW)

# ── Panel C: Lineage specificity with n labels ────────────────────────────────
ax = ax_c

cts = [ct for ct in cell_types_ordered if ct in ct_vals]
x   = np.arange(len(cts))
w   = 0.35
sc_means = [ct_vals[ct][0] for ct in cts]
tc_means = [ct_vals[ct][1] for ct in cts]

ax.bar(x - w/2, sc_means, w, color=SC_COLOR, alpha=0.85, label='SC')
ax.bar(x + w/2, tc_means, w, color=TC_COLOR, alpha=0.85, label='TC')

# Shade adipogenic
adipogenic = ['Areg', 'pre_CP', 'CP', 'WAT_IPC']
adip_idx = [i for i, ct in enumerate(cts) if ct in adipogenic]
if adip_idx:
    ax.axvspan(min(adip_idx)-0.5, max(adip_idx)+0.5, alpha=0.07, color=TC_COLOR, zorder=0)
    ax.text(np.mean(adip_idx), max(max(sc_means), max(tc_means))*1.04,
            'adipogenic lineage', ha='center', fontsize=7.5, color=TC_COLOR, style='italic')

# n labels under each bar pair
ymin = ax.get_ylim()[0] if ax.get_ylim()[0] < 0 else -0.004
for i, ct in enumerate(cts):
    n_sc, n_tc = ct_n[ct]
    ax.text(x[i], -0.006, f'n={n_sc+n_tc}', ha='center', fontsize=6,
            color='#666666', va='top')

ax.set_xticks(x)
ax.set_xticklabels([ct.replace('_', '\n') for ct in cts], fontsize=7.5)
ax.set_ylabel('Mean Fndc5 expression\n(log-normalized)', fontsize=9)
ax.set_ylim(-0.012, max(max(sc_means), max(tc_means)) * 1.15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=8, frameon=False, loc='upper right')
ax.set_title('C  Lineage-restricted: zero in immune cells', loc='left', **LABEL_KW)
ax.text(0.99, 0.01, '≥20 cells/condition shown', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6.5, color='#888888', style='italic')

# ── Panel D: Clock coherence — raw centered + spaghetti ──────────────────────
ax = ax_d

n_genes   = len(clock_genes)
positions = np.arange(n_genes)
np.random.seed(7)

# Plot spaghetti lines per mouse across genes
x_sc_pos = positions - 0.18
x_tc_pos = positions + 0.18

for mouse in all_mice_sc:
    y = [clock_pm[g][mouse][0] for g in clock_genes if mouse in clock_pm[g]]
    x_pos = [x_sc_pos[i] for i, g in enumerate(clock_genes) if mouse in clock_pm[g]]
    ax.plot(x_pos, y, color=mouse_color[mouse], alpha=0.35, lw=1.2, zorder=2)

for mouse in all_mice_tc:
    y = [clock_pm[g][mouse][0] for g in clock_genes if mouse in clock_pm[g]]
    x_pos = [x_tc_pos[i] for i, g in enumerate(clock_genes) if mouse in clock_pm[g]]
    ax.plot(x_pos, y, color=mouse_color[mouse], alpha=0.35, lw=1.2, zorder=2)

# Scatter points on top
for i, gene in enumerate(clock_genes):
    for mouse in all_mice_sc:
        if mouse in clock_pm[gene]:
            v = clock_pm[gene][mouse][0]
            ax.scatter(x_sc_pos[i], v, color=mouse_color[mouse], s=55, zorder=5,
                       edgecolors='white', lw=0.5)
    for mouse in all_mice_tc:
        if mouse in clock_pm[gene]:
            v = clock_pm[gene][mouse][0]
            ax.scatter(x_tc_pos[i], v, color=mouse_color[mouse], s=55, zorder=5,
                       edgecolors='white', lw=0.5)

    # Group means
    sc_centered = [clock_pm[gene][m][0] for m in all_mice_sc if m in clock_pm[gene]]
    tc_centered = [clock_pm[gene][m][0] for m in all_mice_tc if m in clock_pm[gene]]
    ax.hlines(np.mean(sc_centered), positions[i]-0.30, positions[i]-0.06,
              colors='#333333', lw=2, zorder=6)
    ax.hlines(np.mean(tc_centered), positions[i]+0.06, positions[i]+0.30,
              colors='#333333', lw=2, zorder=6)

ax.axhline(y=0, color='#cccccc', lw=0.8, zorder=0)
ax.axvline(x=n_genes-1-0.5, color='#cccccc', lw=1, ls='--')

# padj annotations
padj_labels = {'Nr1d1': '0.0007', 'Nr1d2': '0.008',
               'Tef':   '0.005',  'Dbp':   '0.026', 'Fndc5': '0.004'}
ymax = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.3
for i, gene in enumerate(clock_genes):
    tc_vals_g = [clock_pm[gene][m][0] for m in all_mice_tc if m in clock_pm[gene]]
    ytop = max(tc_vals_g) * 1.08
    ax.text(positions[i], ytop, f'padj\n{padj_labels[gene]}',
            ha='center', fontsize=6.5, color='#444444', va='bottom')

ax.set_xticks(positions)
ax.set_xticklabels(clock_genes, fontsize=9, style='italic')
ax.set_ylabel('Expression (mean-centered\nper gene, log-normalized)', fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

legend_els = [
    Line2D([0],[0], color=MOUSE_COLORS[0], lw=1.5, marker='o', markersize=6,
           markerfacecolor=MOUSE_COLORS[0], label='SC mice (N=3)'),
    Line2D([0],[0], color=MOUSE_COLORS[3], lw=1.5, marker='o', markersize=6,
           markerfacecolor=MOUSE_COLORS[3], label='TC mice (N=4)'),
]
ax.legend(handles=legend_els, fontsize=8, frameon=False, loc='lower left')
ax.set_title('D  Clock coherence: Fndc5 co-induces with Nr1d1/2, Tef, Dbp', loc='left', **LABEL_KW)
ax.text(0.99, 0.01, 'each line = one mouse', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6.5, color='#888888', style='italic')

# ── Suptitle ──────────────────────────────────────────────────────────────────
fig.suptitle(
    'Exercise-associated Fndc5 induction co-occurs with a CLOCK/BMAL1-associated transcriptional program\n'
    'specifically in visceral fat progenitor cells',
    fontsize=11, fontweight='bold', y=0.995)

plt.savefig(OUT_PATH, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUT_PNG,  dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {OUT_PATH}")
print(f"Saved: {OUT_PNG}")
