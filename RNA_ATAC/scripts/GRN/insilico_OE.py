######################
## Import libraries ##
######################

import os
from re import search
from dfply import *
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as s
import math
import seaborn as sns
import importlib
import argparse

import celloracle as co
co.__version__

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

####################
## Load functions ##
####################

def load_adata(adata_file, metadata_file = None, normalise = False, cells = None, cell_column = "cell", features = None, filter_lowly_expressed_genes = False, set_colors = False, keep_counts=False):
    adata = sc.read(adata_file)
    # Convert to sparse matrices
    if not s.sparse.issparse(adata.X):
        adata.X = s.sparse.csr_matrix(adata.X)
    if len(adata.layers.keys())>0:
        for i in list(adata.layers.keys()):
            if not s.sparse.issparse(adata.layers[i]):
                adata.layers[i] = csr_matrix(adata.layers[i])
    if cells is not None:
        tmp = np.mean(np.isin(cells,adata.obs.index.values)==False)
        if tmp<1: print("%.2f%% of cells provided are not observed in the adata, taking the intersect..." % (100*tmp))
        cells = np.intersect1d(cells,adata.obs.index.values)
        adata = adata[cells,:]
    if features is not None:
        adata = adata[:,features]
    if metadata_file is not None:
        metadata = pd.read_table(metadata_file, delimiter="\t", header=0).set_index(cell_column, drop=False)
        metadata = metadata.loc[cells]
        assert np.all(adata.obs.index.isin(metadata[cell_column]))
        # assert np.all(metadata.cell.isin(adata.obs.index))
        assert metadata.shape[0] == adata.shape[0]
        adata.obs = metadata#.reindex(adata.obs.index)
    if filter_lowly_expressed_genes:
        sc.pp.filter_genes(adata, min_counts=10)
    if keep_counts:
        adata.layers["raw"] = adata.X.copy()
    if normalise:
        sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=False)
        sc.pp.log1p(adata)
    if set_colors:
        colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['celltype']))]
        adata.uns['celltype_colors'] = colPalette_celltypes
        colPalette_stages = [opts["stage_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['stage']))]
        adata.uns['stage_colors'] = colPalette_stages
    return adata

###########################
## Load default settings ##
###########################

# if search("BI2404M", os.uname()[1]):
#     exec(open('/Users/argelagr/gastrulation_multiome_10x/settings.py').read())
#     exec(open('/Users/argelagr/gastrulation_multiome_10x/utils.py').read())
# elif search("pebble|headstone", os.uname()[1]):
#    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/settings.py').read())
#    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/utils.py').read())
# else:
#     exit("Computer not recognised")

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--max_expr',              type=str,                required=True,           help='Maximal expression levels')
p.add_argument( '--tf2gene',               type=str,                required=True,           help='TF to gene mapping')
p.add_argument( '--marker_TFs_dir',               type=str,                required=True,           help='Marker TF')
p.add_argument( '--min_chip_score',               type=float,                required=True,           help='Minimal in silico chip score')
p.add_argument( '--max_distance',               type=int,                required=True,           help='Maximal distance')
p.add_argument( '--umap_coord',               type=str,                required=True,           help='UMAP coordinates_MOFA')
p.add_argument( '--umap_coord_RNA',               type=str,                required=True,           help='UMAP coordinates ATAC')
p.add_argument( '--umap_coord_ATAC',               type=str,                required=True,           help='UMAP coordinates RNA')
p.add_argument( '--genes',               type=str,                required=True,   nargs="+",        help='Genes to knock out')
p.add_argument( '--incl_samples',            type=str,    required=True,       default="all",             help='Samples')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
args = p.parse_args()

## START TEST ##
# args = {}
# args["anndata"] = "/data/louisc/Project_Babraham/RNA_ATAC/metacells/anndata_metacells_nodiff.h5ad"
# args["metadata"] = "/data/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
# args["outdir"] = "/data/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/insilico_OE/"
# args["incl_samples"] = "nodiff"
# args["min_chip_score"] = 0.15
# args["max_distance"] = 50000
# args["marker_TFs_dir"] = "/data/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/Marker_TF/"
# args['tf2gene'] = "/data/louisc/Project_Babraham/RNA_ATAC/gene_regulatory_network/metacells/CISBP/TF2gene_after_virtual_chip.txt.gz"
# args['umap_coord'] = "/data/louisc/Project_Babraham/RNA_ATAC/mofa/pdf/umap_" + args["incl_samples"] + ".txt.gz"
# args['umap_coord_RNA'] = "/data/louisc/Project_Babraham/RNA/dimensionality_reduction/batch_correction_None/umap_features2500_pcs15_neigh50_dist0.25_" + args["incl_samples"] + ".txt.gz"
# args['umap_coord_ATAC'] = "/data/louisc/Project_Babraham/ATAC/archR/dimensionality_reduction/umap_PeakMatrix_nfeatures15000_ndims25_" + args["incl_samples"] + ".txt.gz"
# args["genes"] = ["KLF5", "TFCP2L1", "NR5A2", "HNF4G", "ESRRB", "ID3"]
## END TEST ##

opts = {}
opts["celltype_colors"] = ["#DAF7A6","#FFC300","#DDAA05","#C70039","#900C3F","#581845"]

#####################
## Parse arguments ##
#####################

args = vars(args)
print(args)

if not os.path.isdir(args["outdir"]):
    args["outdir"] = Path(args["outdir"])

if args["incl_samples"]=="all":
    args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L"]
elif args["incl_samples"]=="nodiff":
    args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18"]

###################
## Read metadata ##
###################

metadata = (pd.read_table(args["metadata"]) >>
    mask(X["batch"].isin(args["samples"]))
).set_index("sample", drop=False)
print(metadata.shape)
print(metadata.head)

##################
## Load AnnData ##
##################

adata = load_adata(
    adata_file = args['anndata'],
    normalise = True,
    keep_counts = True,
	set_colors = False
)
print(adata) # Check if adata is not yet logtransformed (because it is when loaded as oracle object)
print(adata.X[:5,:5].todense())

# add clusterannotation from metadata
adata.obs["cluster"] = (metadata.loc[adata.obs.index.values] >> select("cluster"))
adata.obs["cluster"] = adata.obs["cluster"].astype('category')

# determin color palette for clusters
colPalette_celltypes = opts["celltype_colors"][:len(np.unique(adata.obs['cluster']))]
adata.uns['cluster'] = colPalette_celltypes

print(adata)

##########################
## Load Expression data ##
##########################

max_expr = pd.read_csv(args["max_expr"], sep="\t")
print(max_expr)
print(max_expr.head())

#######################
## Feature selection ##
#######################

# sc.pp.highly_variable_genes(adata, n_top_genes=25000)
#
# adata = adata[:,adata.var["highly_variable"]]
print(adata)

print(np.isin(np.array(["SOX2"]),adata.var_names))

#########
## PCA ##
#########

sc.tl.pca(adata, n_comps=15, svd_solver='arpack')

######################
## Batch correction ##
######################

# sc.external.pp.harmony_integrate(adata, "stage", basis='X_pca', adjusted_basis='X_pca_harmony')

################
## k-NN graph ##
################

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=15, use_rep="X_pca")

###########################
## Force-directed layout ##
###########################

adata.layers["raw_count"] = adata.layers["raw"]
adata.X = sc.pp.log1p(adata.layers["raw"]).copy()
adata.X[1:5, 1:5]

adata_before_loop = adata

for i in range(1,4):
    print(["MOFA", "RNA", "ATAC"][i-1] + " started!")
    adata = adata_before_loop

    sc.tl.draw_graph(adata, layout='fa', init_pos=None)

    if i == 1:
        adata.obsm["X_draw_graph_fa"] = (pd.read_csv(args["umap_coord"], sep=",", index_col=0).loc[adata.obs.index.values] >> select([0,1])).values
        args["suffix"] = ".pdf"
    elif i == 2:
        adata.obsm["X_draw_graph_fa"] = (pd.read_csv(args["umap_coord_RNA"], sep=",", index_col=2).loc[adata.obs.index.values] >> select([0,1])).values
        args["suffix"] = "_RNA.pdf"
    else:
        adata.obsm["X_draw_graph_fa"] = (pd.read_csv(args["umap_coord_ATAC"], sep=",", index_col=0).loc[adata.obs.index.values] >> select([0, 1])).values
        args["suffix"] = "_ATAC.pdf"

    sc.pl.draw_graph(adata, size=200, color=["cluster"], legend_loc="on data",
                     return_fig=True, palette=colPalette_celltypes).savefig(args["outdir"] + "umap" + args["suffix"])

    #################################
    ## Prepare data for celloracle ##
    #################################

    print("Prepare data for celloracle")

    TFs = []
    for i in range(0,5):
        TF_markers_pandas = pd.read_csv(args["marker_TFs_dir"] + "Marker_TF_cluster" + str(i+1) +".txt", sep="\t") >> (
            mutate(gene=X["Gene"].str.upper())
        )
        #print(TF_markers_pandas.head)
        TFs.append(np.unique(TF_markers_pandas.gene.values))

    TFs = np.unique(np.hstack(TFs))

    tf2gene_df = pd.read_csv(args['tf2gene'],header=0, sep="\t") >> (
        mask(X["chip_score"]>=args["min_chip_score"], X["dist"]<=args["max_distance"]) >>
        mask(X["gene"].str.upper().isin(TFs),X["tf"].str.upper().isin(TFs)) >>
        select([1,3]) >>
        mutate(gene=X["gene"].str.upper())
    )

    print(TFs)
    print(tf2gene_df.shape)
    print(tf2gene_df.head())

    all_genes = np.unique(np.concatenate((tf2gene_df["tf"].values,tf2gene_df["gene"].values)))
    print(len(all_genes))
    # adata.var.index.str.upper().isin()

    adata.var.index = adata.var.index.str.upper()
    adata.var["gene"] = adata.var.index
    adata = adata[:,all_genes]

    tmp = adata.var["gene"].str.upper().isin(tf2gene_df["tf"].values).values
    adata.var["gene"][tmp] = adata.var["gene"][tmp].str.upper()

    adata.var["gene"].str.upper().isin(tf2gene_df["tf"].values).sum()

    tf2gene_df = tf2gene_df[tf2gene_df["tf"].isin(adata.var_names)]
    tf2gene_df.shape

    tf2gene_df = tf2gene_df[tf2gene_df["gene"].isin(adata.var_names)]
    tf2gene_df.shape

    tf2gene_dic = tf2gene_df.groupby('gene')['tf'].apply(lambda x: x.tolist()).to_dict()
    print(len(tf2gene_dic.keys()))
    print(tf2gene_dic.keys())

    ############################
    ## Initiate oracle object ##
    ############################

    oracle = co.Oracle()

    ##################################################
    ## Load gene expression data into oracle object ##
    ##################################################

    adata.var["symbol"] = adata.var.index.values
    #adata.var["isin_top1000_var_mean_genes"] = True

    # create dummy column
    adata.obs["foo"] = "foo"
    adata.obs["foo"] = adata.obs["foo"].astype('category')

    adata.uns["foo_colors"] = np.array(["#9e6762"])

    # check adata object
    adata.obsm

    print(adata)
    print(adata.obsm["X_draw_graph_fa"].shape)
    print(adata.obs["foo"])

    print("Load counts")

    oracle.import_anndata_as_normalized_count(adata=adata, cluster_column_name="foo",
                                              embedding_name="X_draw_graph_fa")

    ####################################
    ## Load TFinfo into oracle object ##
    ####################################

    print("Load TF")

    oracle.import_TF_data(TFdict=tf2gene_dic)

    #oracle.adata.var.loc[:,["isin_TFdict_targets","isin_TFdict_regulators"]]
    print(oracle.adata.var.loc[:,"isin_TFdict_targets"].sum())
    print(oracle.adata.var.loc[:,"isin_TFdict_regulators"].sum())

    #########
    ## PCA ##
    #########

    # Perform PCA
    oracle.perform_PCA()

    # Select important PCs
    fig = plt.figure()
    plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
    fig.savefig(args["outdir"]+"PCA_var_explained.pdf")
    #n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    #plt.axvline(n_comps, c="k")

    ####################
    ## KNN imputation ##
    ####################

    #It creates the attributes:
    # knn (scipy.sparse.csr_matrix): knn contiguity matrix
    # knn_smoothing_w (scipy.sparse.lil_matrix): the weights used for the smoothing

    k = 25
    oracle.knn_imputation(n_pca_dims=15, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)
    oracle.knn.shape
    oracle.adata.layers

    oracle.adata.layers["normalized_count"].todense()[:5,:5]

    #####################
    ## GRN calculation ##
    #####################

    #The next step is constructing a cluster-specific GRN for all clusters.
    #The "get_links" function returns GRNs as a Links object which stores inferred GRNs and the corresponding metadata.
    #You can do network analysis with the Links object.

    links = oracle.get_links(cluster_name_for_GRN_unit="foo", alpha=10, verbose_level=10, test_mode=False)
    links

    links.links_dict["foo"].shape

    links.links_dict["foo"]

    # Exploration
    links.filter_links(p=0.0001, weight="coef_abs")
    links.links_dict["foo"].shape

    # Save
    args['links_outfile'] = args['outdir'] + 'test.celloracle.links'
    args['oracle_outfile'] = args['outdir'] + "test.celloracle.oracle"

    print("GRN")
    print(links)
    print(links.links_dict)
    print(links.links_dict["foo"].shape)
    print(links.links_dict["foo"])

    #############
    ## Explore ##
    #############

    for gene in TFs:
        print(gene)

        # Explore links for a gene
        # print(links.links_dict["foo"] >> mask(X.source==gene) >> arrange("coef_mean"))
        #
        oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
        #
        oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

        if gene not in oracle.active_regulatory_genes:
            print("Skipped due to insufficient connections in GRN!")
            continue
        #
        # sc.pl.draw_graph(oracle.adata, color=[gene], size=200, layer="normalized_count", cmap="viridis",save=True)
        #
        cutoff_2fold = float(sc.pp.log1p(max_expr["max_expr"][max_expr["gene"] == gene] * 2))
        print("cutoff is " + str(cutoff_2fold))
        oracle.simulate_shift(perturb_condition={gene: cutoff_2fold}, n_propagation=3)
        #
        oracle.estimate_transition_prob(n_neighbors=15, knn_random=True, sampled_fraction=1)
        #
        oracle.calculate_embedding_shift(sigma_corr = 0.05)
        #
        oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=15)
        fig=plt.figure()
        plt.hist(oracle.total_p_mass)
        plt.show()
        fig.savefig(args["outdir"] + gene + "_histogram_p_mass" + args["suffix"])
        #
        oracle.suggest_mass_thresholds(n_suggestion=8)
        #
        oracle.calculate_mass_filter(min_mass=0.5, plot=True)
        #
        scale_simulation = 15
        fig= plt.figure()
        oracle.plot_simulation_flow_on_grid(scale=scale_simulation)
        fig.savefig(args['outdir'] + "Insilio_OE_" + gene + "_flow" + args["suffix"])
        #
        cell_colors = [opts["celltype_colors"][i-1] for i in oracle.adata.obs["cluster"].values]
        fig, ax = plt.subplots(figsize=[10, 8])
        ax.scatter(oracle.embedding[:,0], oracle.embedding[:,1], c=cell_colors, alpha=0.45)
        oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
        plt.savefig(args['outdir'] + "Insilio_OE_" + gene + "_flow_col" + args["suffix"])

#############################
## Create Completion token ##
#############################

f = open(args["outdir"]+"completed_"+args["incl_samples"]+".txt", "x")