import scanpy as sc
import pandas as pd
import os  # Import the os module

# Define the results file path
results_file = "write/pbmc3k.h5ad"
pbmc = sc.datasets.pbmc68k_reduced()

# Ensure the 'write' directory exists
if not os.path.exists('write'):
    os.makedirs('write')

#------------------------------------------------------------------------------------------
#Preprocessing and clustering 3k PBMCs (legacy workflow)
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")

results_file = "write/pbmc3k.h5ad"  # the file that will store the analysis results

adata = sc.read_10x_mtx(
    "C:/Users/Alice/Documents/APPPROJECT/tests/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx
adata
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#Preprocessing
sc.pl.highest_expr_genes(adata, n_top=20)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)


# A violin plot of some of the computed quality measures:
# the number of genes expressed in the count matrix
# the total counts per cell
# the percentage of counts in mitochondrial genes
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

#Remove cells that have too many mitochondrial genes expressed or too many total counts
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")

#Actually do the filtering by slicing the AnnData object.
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

#Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.normalize_total(adata, target_sum=1e4)

#Logarithmize the data:
sc.pp.log1p(adata)

#Identify highly-variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)

adata.raw = adata

#Actually do the filtering
adata = adata[:, adata.var.highly_variable]

#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

#Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)
print("-----------------------preprocessing finished")
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#Principal component analysis
#We can make a scatter plot in the PCA coordinates, but we will not use that later on.
print("-----------------------PCA Starting")
sc.tl.pca(adata)
sc.pl.pca(adata, color="CST3")
print("--------------------------scatterplot PCA finished")
#contribution of single PCs to the total variance in the data
sc.pl.pca_variance_ratio(adata, log=True)
print("--------------------------plot variance ratio finished")

#Save the result.
adata.write(results_file)
adata
print("-------------------------results saved, PCA Ending")
#------------------------------------------------------------------------------------------


#---------------------------------
#Computing the neighborhood graph
print("-------------------------computing the neighborhood graph")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
print("-------------------------END computing the neighborhood graph")

#---------------------------------
#Embedding the neighborhood graph
print("-------------------------Embedding the neighborhood graph")
sc.tl.leiden(adata)
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')

sc.tl.umap(adata)

sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"])

sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"], use_raw=False)
print("-------------------------END Embedding the neighborhood graph")
#------------------------------------------------------------------------------------------


#---------------------------------
#Clustering the neighborhood graph
print("-------------------------Clustering the neighborhood graph")
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=None,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)

sc.pl.umap(adata, color=["leiden", "CST3", "NKG7"])
adata.write(results_file)
print("-------------------------END Clustering the neighborhood graph")
#------------------------------------------------------------------------------------------


#---------------------------------
#Finding marker genes
print("-------------------------Finding marker genes")

#ranking for the highly differential genes
sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
print("-------------------------Ranking genes groups finished")

sc.settings.verbosity = 2  # reduce the verbosity
print("-------------------------Verbosity reduction finished")

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
print("-------------------------Ranking genes groups #2 finished")

#Save the result.
adata.write(results_file)
print("-------------------------Result saved")

#rank genes using logistic regression
sc.tl.rank_genes_groups(adata, "leiden", method="logreg", max_iter=1000)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
print("-------------------------Ranking genes groups using logistic regression finished")

#define a list of marker genes for later reference
marker_genes = [
    *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
    *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
    *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]
print("-------------------------definition of markeer genes finished")

#Reload the object that has been save with the Wilcoxon Rank-Sum test result.
adata = sc.read(results_file)

#how the 10 top ranked genes per cluster 0, 1, …, 7 in a dataframe.
pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)

#Get a table with the scores and groups.
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "pvals"]
    }
).head(5)

#Compare to a single cluster:
sc.tl.rank_genes_groups(adata, "leiden", groups=["0"], reference="1", method="wilcoxon")
sc.pl.rank_genes_groups(adata, groups=["0"], n_genes=20)

#f we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin.
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)

#Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups):
adata = sc.read(results_file)

sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)

#If you want to compare a certain gene across groups, use the following.
sc.pl.violin(adata, ["CST3", "NKG7", "PPBP"], groupby="leiden")

#Actually mark the cell types.
new_cluster_names = [
    "CD4 T",
    "B",
    "FCGR3A+ Monocytes",
    "NK",
    "CD8 T",
    "CD14+ Monocytes",
    "Dendritic",
    "Megakaryocytes",
]
adata.rename_categories("leiden", new_cluster_names)

sc.pl.umap(
    adata, color="leiden", legend_loc="on data", title="", frameon=False, save=".pdf"
)

#Now that we annotated the cell types, let us visualize the marker genes.
sc.pl.dotplot(adata, marker_genes, groupby="leiden");

#There is also a very compact violin plot.
sc.pl.stacked_violin(adata, marker_genes, groupby="leiden");

adata
# `compression='gzip'` saves disk space, and slightly slows down writing and subsequent reading
adata.write(results_file, compression="gzip")
adata.raw.to_adata().write("./write/pbmc3k_withoutX.h5ad")

print("-------------------------END Finding marker genes")
#------------------------------------------------------------------------------------------

