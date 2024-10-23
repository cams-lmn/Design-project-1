######################
## Import libraries ##
######################

import os
from re import search
import scvelo as scv

import io
import numpy as np
import anndata as anndata
import scipy as s
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import csv
from dfply import *
# from IPython.display import display
from pathlib import Path
import argparse

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
                adata.layers[i] = s.sparse.csr_matrix(adata.layers[i])
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

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata_scvelo',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--incl_samples',            type=str,    required=True,       default="all",             help='Samples')
p.add_argument( '--umap_coord',               type=str,                required=True,           help='UMAP coordinates MOFA')
p.add_argument( '--umap_coord_RNA',               type=str,                required=True,           help='UMAP coordinates RNA')
p.add_argument( '--umap_coord_ATAC',               type=str,                required=True,           help='UMAP coordinates ATAC')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
args = p.parse_args()

# ## START TEST ##
# args = {}
# args["anndata_scvelo"] = "/data/homes/louisc/Project_Babraham/RNA/velocyto/anndata_scvelo.h5ad"
# args["incl_samples"] = "nodiff"
# args['umap_coord'] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/pdf/umap_" + args["incl_samples"] + ".txt.gz"
# args['umap_coord_RNA'] = "/data/homes/louisc/Project_Babraham/RNA/dimensionality_reduction/batch_correction_None/umap_features2500_pcs15_neigh50_dist0.25_" + args["incl_samples"] + ".txt.gz"
# args['umap_coord_ATAC'] = "/data/homes/louisc/Project_Babraham/ATAC/archR/dimensionality_reduction/umap_PeakMatrix_nfeatures15000_ndims25_" + args["incl_samples"] + ".txt.gz"
# args["metadata"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
# args["outdir"] = "/data/homes/louisc/Project_Babraham/RNA/velocyto/"
# ## END TEST ##

# io = {}
# io["cellranger_output"] = "/data/homes/louisc/Project_Babraham/results_sierra"

# convert args to dictionary
args = vars(args)

opts = {}
opts["celltype_colors"] = ["#DAF7A6","#FFC300","#DDAA05","#C70039","#900C3F","#581845"]

####################
## Define options ##
####################

sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(8, 7), facecolor='white')
sc.settings.figdir = args["outdir"]

scv.set_figure_params('scvelo')  # for beautified visualization

########################
## Load cell metadata ##
########################

metadata = (pd.read_table(args["metadata"],index_col=0))
print(metadata.shape)
print(metadata.head)

metadata.index.value_counts()
metadata["cluster"].value_counts()

#########################
## Load anndata object ##
#########################

adata = load_adata(
    adata_file=args["anndata_scvelo"],
    cells=metadata.index.values,
    metadata_file=args["metadata"],
    cell_column = "sample",
    normalise=True,
    filter_lowly_expressed_genes=False,
)
adata

# add clusterannotation from metadata
adata.obs["cluster"] = (metadata.loc[adata.obs.index.values] >> select("cluster"))
adata.obs["cluster"] = adata.obs["cluster"].astype('category')

# determin color palette for clusters
colPalette_celltypes = opts["celltype_colors"][:len(np.unique(adata.obs['cluster']))]
adata.uns['cluster'] = colPalette_celltypes

adata.obs["cluster"].value_counts()

print(adata)

###########################
## Force-directed layout ##
###########################

sc.tl.draw_graph(adata, layout='fa', init_pos=None)

# RNA umap
adata.obsm["X_draw_graph_fa"] = (pd.read_csv(args["umap_coord_RNA"], sep=",", index_col=2).loc[adata.obs.index.values] >> select([0,1])).values
adata.obsm["X_draw_graph_fa_RNA"] = (pd.read_csv(args["umap_coord_RNA"], sep=",", index_col=2).loc[adata.obs.index.values] >> select([0,1])).values
sc.pl.draw_graph(adata, size=200, color=["cluster"], legend_loc="on data",
                 return_fig=True, palette=colPalette_celltypes).savefig(args["outdir"] + "umap_RNA.pdf")

# ATAC umap
adata.obsm["X_draw_graph_fa"] = (pd.read_csv(args["umap_coord_ATAC"], sep=",", index_col=0).loc[adata.obs.index.values] >> select([0,1])).values
adata.obsm["X_draw_graph_fa_ATAC"] = (pd.read_csv(args["umap_coord_ATAC"], sep=",", index_col=0).loc[adata.obs.index.values] >> select([0,1])).values
sc.pl.draw_graph(adata, size=200, color=["cluster"], legend_loc="on data",
                 return_fig=True, palette=colPalette_celltypes).savefig(args["outdir"] + "umap_ATAC.pdf")

# Mofa umap
adata.obsm["X_draw_graph_fa"] = (pd.read_csv(args["umap_coord"], sep=",", index_col=0).loc[adata.obs.index.values] >> select([0,1])).values
sc.pl.draw_graph(adata, size=200, color=["cluster"], legend_loc="on data",
                 return_fig=True, palette=colPalette_celltypes).savefig(args["outdir"] + "umap.pdf")

################
## k-NN graph ##
################

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=25)

################
## k-NN graph ##
################

scv.tl.velocity_graph(adata, approx=True)

for vkey in ["velocity"]:
    adata.layers[vkey] = adata.layers[vkey].todense()

# RNA
scv.pl.velocity_embedding_grid(
    adata, basis='X_draw_graph_fa_RNA', color=['cluster'], size=20, alpha=0.45,
    arrow_length=5, arrow_size=2, dpi=500, legend_loc="none", save=args["outdir"] + "layout_velocity_RNA.pdf", figsize=(5, 5)
)
scv.pl.velocity_embedding_stream(adata, basis='X_draw_graph_fa_RNA', color=["cluster"],alpha=0.45,
    legend_loc='right margin', save=args["outdir"] + "layout_velocity_stream_RNA.png", linewidth=1, legend_fontsize=8, dpi=500)

# ATAC
scv.pl.velocity_embedding_grid(
    adata, basis='X_draw_graph_fa_ATAC', color=['cluster'], size=20, alpha=0.45,
    arrow_length=5, arrow_size=2, dpi=500, legend_loc="none", save=args["outdir"] + "layout_velocity_ATAC.pdf", figsize=(5, 5)
)
scv.pl.velocity_embedding_stream(adata, basis='X_draw_graph_fa_ATAC', color=["cluster"],
    legend_loc='right margin', save=args["outdir"] + "layout_velocity_stream_ATAC.png", linewidth=1, legend_fontsize=8, dpi=500)

# MOFA
scv.pl.velocity_embedding_grid(
    adata, basis='X_draw_graph_fa', color=['cluster'], size=20, alpha=0.45,
    arrow_length=5, arrow_size=2, dpi=500, legend_loc="none", save=args["outdir"] + "layout_velocity.pdf", figsize=(5, 5)
)
scv.pl.velocity_embedding_stream(adata, basis='X_draw_graph_fa', color=["cluster"],
    legend_loc='right margin', save=args["outdir"] + "layout_velocity_stream.png", linewidth=1, legend_fontsize=8, dpi=500)

#############################
## Create Completion token ##
#############################

f = open(args["outdir"]+"scvelo_analysis_cells_completed_"+args["incl_samples"]+".txt", "x")