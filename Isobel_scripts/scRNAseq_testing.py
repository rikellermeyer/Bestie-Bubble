#!/Users/pfb2024/mamba/envs/bestie_bubble/bin/python

import scanpy as sc
import anndata as ad
import pooch


sc.settings.set_figure_params(dpi=50, facecolor="white")

EXAMPLE_DATA = pooch.create(
    path=pooch.os_cache("scverse_tutorials"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
)
EXAMPLE_DATA.load_registry_from_doi()

samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}

for sample_id, filename in samples.items():
    path = EXAMPLE_DATA.fetch(filename)
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample") #adata is matrix with slots for extra things
adata.obs_names_make_unique()
## print(adata.obs["sample"].value_counts())
## adata
print(f'before anything: {adata}')
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

#plot number of genes in matrix, counts per cell, % mitochondrial
## sc.pl.violin(
##     adata,
##     ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
##     jitter=0.4,
##     multi_panel=True,
## )

##joint QC plot (total counts on x, counts per cell on y, % mitochondrial as color)
#sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.scrublet(adata, batch_key="sample")

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample") #project, expects logarithmized data except when flavor ='seurat_v3'

##sc.pl.highly_variable_genes(adata) #plot highly variable genes

print(f'before pca: {adata}')
sc.tl.pca(adata)
print(f'after pca, before neighborhood: {adata}')

##sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True) #plot contribution of single PCs to total variance

sc.pp.neighbors(adata) #compute neighborhood graph

print(f'after neighborhood, before umap: {adata}')

sc.tl.umap(adata) #embed in 2 dimensions for visualization

print(f'after umap, before clustering: {adata}')

sc.pl.umap( #plotting UMAP, colored by sample
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

print(f'after clustering: {adata}')

sc.pl.umap(adata, color=["leiden"])
