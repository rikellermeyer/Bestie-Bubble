#!/Users/pfb2024/mamba/envs/bestie_bubble/bin/python

import scanpy as sc
import anndata as ad

sc.settings.set_figure_params(dpi=250, facecolor="white")

with open("1000_data_trimmed2.csv") as your_data:
    adata = ad.read_csv(your_data, delimiter='\t')

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)

sc.tl.pca(adata)

sc.pp.neighbors(adata)

sc.tl.umap(adata) 

sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

sc.pl.umap(adata, color=["leiden"], size=40)