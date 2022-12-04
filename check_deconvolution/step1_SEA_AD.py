
import scanpy as sc
import numpy as np
import anndata as ad
import glob

import re

# modified based on github.com/scverse/scanpy/scanpy/preprocessing/_utils.py
def get_mean_var(X, axis=0):
    mean = np.mean(X, axis=axis)
    mean_sq = np.multiply(X, X).mean(axis=axis)
    var = mean_sq - mean**2
    var *= X.shape[axis] / (X.shape[axis] - 1)
    return mean, var

# --------------------------------------------------
# list file names
# --------------------------------------------------

work_dir = '/Users/wsun/research/data/SEA-AD/'
d = '/Users/wsun/research/data/SEA-AD/*.h5ad'

files = glob.glob(d)
len(files)


# --------------------------------------------------
# read in data
# --------------------------------------------------

for f1 in files:
    print(f1)
    cell_type = re.findall("SEA-AD/(\S+).h5ad", f1)[0]
    cell_type

    adata = ad.read_h5ad(f1)
    adata.obs

    list(adata.obs.columns.values)
    adata.obs.donor_id.value_counts()
    adata.obs.cell_type.value_counts()
    adata.obs.disease.value_counts()

    adata.var.feature_is_filtered.value_counts()
    adata.var.feature_biotype.value_counts()
    adata.var.feature_reference.value_counts()
    adata.var.iloc[0:5,1:3]

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False,
                               inplace=True, use_raw=True)

    list(adata.obs.columns.values)
    list(adata.var.columns.values)

    adata.obs.total_counts.describe()
    adata.var.total_counts.describe()
    adata.raw.X[5:8, 3:18].toarray()
    adata.X[5:8, 3:18].toarray()

    edat = np.expm1(adata.X.toarray())
    edat.shape
    edat[5:8, 3:18]

    mean, var = get_mean_var(edat, axis=0)

    df = adata.var
    df['mean_rd_normalized_count'] = mean
    df['var_rd_normalized_count'] = var

    df.to_csv(work_dir+cell_type+'_genes_info.csv')

