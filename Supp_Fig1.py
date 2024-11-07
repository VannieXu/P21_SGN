import scanpy as sc
import pandas as pd
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import copy

vir = copy(mpl.cm.viridis)
vir.set_under("lightgray")

import numpy as np

sc.settings.verbosity = 3            
sc.logging.print_header()
sc.settings.set_figure_params(figsize = [10,8],dpi =120, facecolor='white')

adata = sc.read_h5ad('P21_raw.h5ad')

adata.layers['plot'] = np.log1p(adata.X)
adata.layers['count'] = adata.X.copy()

adata.var['n_cell'] = np.ravel(adata.X.getnnz(axis = 0))
adata.var['total_count'] = np.ravel(adata.X.sum(axis = 0))

adata.obs['n_gene'] = np.ravel(adata.layers['count'].getnnz(axis = 1))
adata.obs['total_UMI'] = np.ravel(adata.layers['count'].sum(axis = 1))

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata)

sc.tl.pca(adata,svd_solver = 'auto')
sc.pl.pca_variance_ratio(adata,n_pcs=50)

sc.pp.neighbors(adata,metric = 'sqeuclidean',n_pcs = 20,n_neighbors = 50)

sc.tl.umap(adata,min_dist = 0.3)

genes = ['Pou4f3','Skp1a','Mog','Mpz','Esam','Vmo1','C1qa','Cd79a','Tubb3','Gria2']
sc.pl.umap(adata,color = genes,layer = 'plot',color_map = vir,vmin = 0.0000001,save = '_marker.pdf',frameon = False,
          ncols = 3)

sc.pl.umap(adata,color = 'total_UMI',vmax = 30000,color_map = 'RdYlGn',save = '_umi.pdf',frameon = False)
sc.pl.umap(adata,color = 'n_gene',color_map = 'plasma',vmax = 10000,save = 'ngene.pdf',frameon = False)
