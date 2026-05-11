import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import anndata as ad

WORKDIR = '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise'
TMP = os.path.join(WORKDIR, 'seurat_tmp')

print('Loading counts matrix...')
X = scipy.io.mmread(os.path.join(TMP, 'counts.mtx'))
X = scipy.sparse.csr_matrix(X).T  # transpose: genes x cells -> cells x genes
print(f'  shape: {X.shape}')

print('Loading gene and barcode names...')
genes = pd.read_csv(os.path.join(TMP, 'genes.txt'), header=None)[0].values
barcodes = pd.read_csv(os.path.join(TMP, 'barcodes.txt'), header=None)[0].values

print('Loading metadata...')
obs = pd.read_csv(os.path.join(TMP, 'metadata.csv'), index_col=0)

var = pd.DataFrame(index=genes)

print('Building AnnData...')
adata = ad.AnnData(X=X, obs=obs, var=var)
adata.obs_names = barcodes
adata.var_names = genes

# Add normalized data as a layer if available
data_mtx_path = os.path.join(TMP, 'data.mtx')
if os.path.exists(data_mtx_path):
    print('Loading normalized data...')
    data = scipy.io.mmread(data_mtx_path)
    data = scipy.sparse.csr_matrix(data).T
    adata.layers['counts'] = adata.X.copy()
    adata.X = data
    print('  X = log-normalized, layers["counts"] = raw counts')

# Add dimensional reductions
for red_name in ['pca', 'tsne', 'umap']:
    csv_path = os.path.join(TMP, f'{red_name}.csv')
    if os.path.exists(csv_path):
        print(f'Loading reduction: {red_name}...')
        emb = pd.read_csv(csv_path, index_col=0).values
        adata.obsm[f'X_{red_name}'] = emb

out_path = os.path.join(WORKDIR, 'GSE183288_Single_cell_atlas.h5ad')
print(f'Writing to {out_path}...')
adata.write_h5ad(out_path)
print('Done!')
print(adata)
