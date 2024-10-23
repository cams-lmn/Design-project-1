######################
## Import libraries ##
######################

import os
from re import search
import SEACells

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
import random
import argparse
import sys
import palantir

####################
## Load functions ##
####################

def load_adata(adata_file, metadata_file = None, normalise = False, cells = None, cell_column = "cell", features = None, filter_lowly_expressed_genes = False, set_colors = False, keep_counts=False):
    adata = sc.read(adata_file)
    # Convert to sparse matrices
    if not s.sparse.issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    if len(adata.layers.keys())>0:
        for i in list(adata.layers.keys()):
            if not issparse(adata.layers[i]):
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
#    exec(open('/Users/argelagr/gastrulation_multiome_10x/settings.py').read())
#     exec(open('/Users/argelagr/gastrulation_multiome_10x/utils.py').read())
# elif search("pebble|headstone", os.uname()[1]):
#     exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/settings.py').read())
#     exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/utils.py').read())
# else:
#     exit("Computer not recognised")

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='')
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--samples',            type=str, nargs="+",              required=True,            help='samples to use')
p.add_argument( '--percent_metacells',            type=float,              default=0.05,             help='Number of metacells (as a fraction of the total number of cells)')
p.add_argument( '--n_features',            type=int,              default=1500,             help='Number of features')
p.add_argument( '--n_pcs',            type=int,              default=15,             help='Number of PCs')
p.add_argument( '--seed',                  type=int,                default=42,               help='Random seed')
p.add_argument( '--seed2',                  type=int,                default=4242,               help='Random seed 2')
p.add_argument( '--seed3',                  type=int,                default=424242,               help='Random seed 3')
p.add_argument( '--n_iter',       type=int,              default=50,              help='Number of iterations')
args = p.parse_args()

# ## START TEST ##
# args = {}
# args["anndata"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/anndata_nodiff.h5ad"
# args["metadata"] =  "/data/homes/louisc/Project_Babraham/RNA/mapping/sample_metadata_after_doublets.txt.gz"
# args["outdir"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/"
# args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18"]
# args["percent_metacells"] = 0.0133
# args["n_features"] = 1500
# args["n_pcs"] = 15
# args["seed"] = 42
# args["seed2"] = 4242
# args["seed3"] = 424242
# args["n_iter"] = 50
# ## END TEST ##

# convert args to dictionary
args = vars(args)
print(args)

# set seed
random.seed(args["seed"])

#####################
## Parse arguments ##
#####################

# I/O
# io["pca_rna"] = io["basedir"] + "/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz"
# io["pca_atac"] = io["basedir"] + "/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh45_dist0.45.txt.gz"

if not os.path.isdir(args["outdir"]):
    os.makedirs(args["outdir"])
    args["outdir"] = Path(args["outdir"])


sc.settings.figdir = args["outdir"] + "pdf"

if not os.path.isdir(sc.settings.figdir):
    os.makedirs(sc.settings.figdir)

# if isinstance(args["samples"],list):
#   if args["samples"][0]=="all":
#     args["samples"] = opts["samples"]
#   else:
#     assert set(args["samples"]).issubset(opts["samples"])

# else:
#   print('args["samples"] has to be a list')


###################
## Load metadata ##
###################

metadata = (pd.read_table(args["metadata"]) >>
	# mask(X["pass_rnaQC"]==True, X["pass_atacQC"]==True, X["doublet_call"]==False, X["celltype"].isin(opts["celltypes"])) >>
	# mask(X["pass_rnaQC"]==True, X["doublet_call"]==False, X["celltype"].isin(opts["celltypes"])) >>
	mask(X["pass_rnaQC"] == True, X["doublet_call"] == False) >>
	mask(X["sample"].isin(args["samples"]))
).set_index("cell", drop=False)

# TO-DO: USE ONLY WT CELLS 

print(metadata.shape)
print(metadata.head())

##################
## Load AnnData ##
##################

adata = load_adata(
	adata_file = args["anndata"], 
	metadata_file = args["metadata"], 
	cells = metadata.index.values, 
	normalise = True, 
	keep_counts = True,
	filter_lowly_expressed_genes = True, 
	set_colors = False
)

# Set colors
# adata.obs = adata.obs.rename(columns={"celltype.mapped":"celltype"})
# colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['celltype']))]
# adata.uns['celltype_colors'] = colPalette_celltypes
# colPalette_stages = [opts["stages_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
# adata.uns['stage_colors'] = colPalette_stages


#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=args["n_features"])

##############################
## Dimensionality reduction ##
##############################

random.seed(args["seed3"])

# Load precomputed PCA coordinates
# pca_mtx = pd.read_csv(io["pca_rna"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
# pca_mtx = pd.read_csv(io["pca_atac"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
# adata.obsm["X_pca"] = pca_mtx

# Run PCA
sc.tl.pca(adata, svd_solver='arpack', n_comps=args["n_pcs"])

# Plot PCA
# sc.pl.pca(adata, components=[1,2], color=["celltype","stage"], size=25, legend_loc=None)

# Build kNN graph
sc.pp.neighbors(adata, n_neighbors=25, use_rep='X_pca')

# Run UMAP
sc.tl.umap(adata, min_dist=0.5, n_components=2)

# Plot UMAP
sc.pl.umap(adata, size=25, legend_loc=None, save="_umap_cells.pdf")

########################
## Fit metacell model ##
########################

random.seed(args["seed3"])

n_metacells = round(args["percent_metacells"] * adata.shape[0])

print("Fitting SEACells with %d metacells..." % (n_metacells))

model = SEACells.core.SEACells(adata,
				  build_kernel_on = 'X_pca',
				  n_SEACells = n_metacells,
				  n_waypoint_eigs = 10,
				  convergence_epsilon = 1e-4)

model.construct_kernel_matrix()

model.initialize_archetypes()

SEACells.plot.plot_initialization(adata, model, save_as=args["outdir"] + "pdf/initialisation.pdf")

model.fit(max_iter=args["n_iter"])

print(pd.value_counts(adata.obs[["SEACell"]].values.transpose().tolist()[0]))

adata.obs[["SEACell_old"]] = adata.obs[["SEACell"]]

mask_SEACells = np.isin(adata.obs[['SEACell']].index.values,model.get_hard_archetypes().values)
adata.obs[['SEACell']].values[mask_SEACells]

for i in range(0,315):
    SEACell_i = adata.obs[["SEACell"]].values[mask_SEACells][i][0]
    mask_SEACell_i = adata.obs[["SEACell"]].values==SEACell_i
    mask_SEACell_i = mask_SEACell_i.transpose()
    mask_SEACell_i = mask_SEACell_i.tolist()[0]
    SEACell_i_new = adata.obs[["SEACell"]].index.values[mask_SEACell_i][np.isin(adata.obs[["SEACell"]].index.values[mask_SEACell_i],
                                                                                                                model.get_hard_archetypes().values)][0]
    adata.obs.loc[adata.obs[["SEACell"]].index.values[mask_SEACell_i],"SEACell"] = SEACell_i_new

print(pd.value_counts(adata.obs[["SEACell"]].values.transpose().tolist()[0]))

#######################
## Plot model output ##
#######################

model.plot_convergence(save_as=args["outdir"] + "pdf/model_convergence.pdf")

SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=False, save_as=args["outdir"] + "pdf/umap_highlight_metacells.pdf")

################################################################
## Aggregate counts and plot trajectory at the metacell level ##
################################################################

adata_metacells = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', celltype_label="cell", summarize_layer='raw')
adata_metacells.uns = adata.uns
adata_metacells.obs = (adata.obs.loc[adata_metacells.obs.cell] >>
    select(["cell","sample"])
)
# adata_metacells.obs = (adata.obs.loc[adata_metacells.obs.index] >>
#     select(["cell","sample","cluster"])
# )
sc.pp.normalize_total(adata_metacells)
sc.pp.log1p(adata_metacells)
sc.pp.highly_variable_genes(adata_metacells, n_top_genes=2500)
sc.tl.pca(adata_metacells, n_comps=25)
sc.pp.neighbors(adata_metacells, n_neighbors=25, use_rep='X_pca')
sc.tl.umap(adata, min_dist=0.5, n_components=2)
sc.pl.umap(adata, size=25, legend_loc=None, save="_umap_metacells.pdf")

##########
## Save ##
##########

to_save = adata.obs[['SEACell','SEACell_old']].reset_index()
to_save.columns = ["cell","metacell","metacell_alt"]

outfile = args["outdir"] + "cell2metacell_assignment_nodiff.txt.gz"
to_save.to_csv(outfile, sep="\t", header=True, index=False)

adata_metacells.write_h5ad(args["outdir"] + "anndata_metacells_nodiff.h5ad")
