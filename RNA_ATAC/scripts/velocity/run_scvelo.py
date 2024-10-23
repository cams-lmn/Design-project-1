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
        adata.X = csr_matrix(adata.X)
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
#     exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/settings.py').read())
#     exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/utils.py').read())
# else:
#     exit("Computer not recognised")

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--incl_samples',            type=str,  required=True,       default="all",             help='Samples')
p.add_argument( '--ncores',            type=int,              default=1,             help='Number of cores')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
args = p.parse_args()

## START TEST ##
# args = {}
# args["anndata"] = "/data/homes/louisc/Project_Babraham/RNA/velocyto/anndata_velocyto.h5ad"
# args["metadata"] = "/data/homes/louisc/Project_Babraham/RNA/mapping/sample_metadata_after_doublets.txt.gz"
# args["incl_samples"] = "all"
# args["ncores"] = 15
# args["outdir"] = "/data/homes/louisc/Project_Babraham/RNA/velocyto"
## END TEST ##

# convert args to dictionary
args = vars(args)

#####################
## Parse arguments ##
#####################

if not os.path.isdir(args["outdir"]): os.makedirs(args["outdir"])
# if not os.path.isdir(os.path.dirname(args["outfile"])):
#     os.makedirs(os.path.dirname(args["outfile"]))

if args["incl_samples"]=="all":
    args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L"]
elif args["incl_samples"]=="nodiff":
    args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18"]

print(args)

###################
## Load metadata ##
###################

print("Loading metadata...")

metadata = (pd.read_table(args["metadata"]) >>
    mask(X.pass_rnaQC==True, X.doublet_call==False, X["sample"].isin(args["samples"]))
).set_index("cell", drop=False)

print(metadata.head())

##################
## Load anndata ##
##################

print("Loading anndata...")

adata = load_adata(
    adata_file = args["anndata"],
    metadata_file = args["metadata"],
    normalise = True,
    cells = metadata.index.values
)

print(adata)

assert "spliced" in list(adata.layers.keys())
assert "unspliced" in list(adata.layers.keys())

############
## scVelo ##
############

# Plot exonic vs intronic read proportions
scv.pl.proportions(adata, save = args["outdir"]+"/proportion_reads.pdf")

# Gene filter
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# Calculate means variances among nearest neighbors in PCA space
scv.pp.moments(adata, n_pcs=50, n_neighbors=25)

# Recover the full splicing kinetics of specified genes.
# For each gene the model infers transcription rates, splicing rates, degradation rates, as well as cell-specific latent time and transcriptional states
scv.tl.recover_dynamics(adata, n_jobs=args["ncores"])

# Estimate velocities per gene
scv.tl.velocity(adata, mode="dynamical")

##########
## Save ##
##########

print("Saving output...")

# Rename index
adata.obs.index.name = "cells"

# delete unnecesary layers
del adata.layers["spliced"]

# cast some layers to float32 to reduce disk
adata.layers["fit_t"] = adata.layers["fit_t"].astype(np.float32)
adata.layers["fit_tau"] = adata.layers["fit_tau"].astype(np.float32)
adata.layers["fit_tau_"] = adata.layers["fit_tau_"].astype(np.float32)
adata.layers["velocity"] = adata.layers["velocity"].astype(np.float32)
adata.layers["velocity_u"] = adata.layers["velocity_u"].astype(np.float32)

# Logical columns need to be stored as str to avoid hdf5 complaining
adata.obs["pass_rnaQC"] = adata.obs["pass_rnaQC"].astype(str)
adata.obs["doublet_call"] = adata.obs["doublet_call"].astype(str)

adata.write_h5ad(args["outdir"]+"/anndata_scvelo_" + args["incl_samples"] + ".h5ad", compression="gzip")

# peek_shown = ["FAM64A", "RNASEH2B", "ELAVL4", "DCX", "STMN2", "GRIA3"]
# plt.figure(None, (7,6.))
# gs = plt.GridSpec(3,2)
# for n, gene in enumerate(peek_shown):
#     i = np.where(vlm.ra["Gene"] == gene)
#     ax = plt.subplot(gs[n])
#     plt.scatter(pc_obj.arclength[pc_obj.ixsort], vlm.Ux_sz[i, pc_obj.ixsort],
#                 alpha=0.7, c=np.array([0,159,193])/255, s=5, label="unspliced")
#     plt.scatter(pc_obj.arclength[pc_obj.ixsort], vlm.Sx_sz[i, pc_obj.ixsort]*vlm.gammas[i],
#                 alpha=0.7, c=np.array([251, 172, 71])/255, s=5, label="spliced")
#     m = 0 #np.minimum(np.min(vlm.Ux_sz[i,:]), np.min(vlm.Sx_sz[i,:]*vlm.gammas[i]))
#     M = np.maximum(np.max(vlm.Ux_sz[i,:]), np.max(vlm.Sx_sz[i,:]*vlm.gammas[i]))
#     plt.ylim(m - 0.07*(M-m), M + 0.07*(M-m))
#     plt.ylabel(gene)
#     plt.yticks([m,0.5*(m+M),M], [f"{m:.2f}", "", f"{M:.2f}"])
#     p = np.min(pc_obj.arclength[pc_obj.ixsort])
#     P = np.max(pc_obj.arclength[pc_obj.ixsort])
#     plt.xticks(np.linspace(p,P,5), [f"{p:.0f}", "","","", f"{P:.0f}"])
#     # Hide the right and top spines
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     # Only show ticks on the left and bottom spines
#     ax.yaxis.set_ticks_position('left')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.spines['left'].set_bounds(m, M)
#     ax.spines['bottom'].set_bounds(p, P)
#     if n == 1:
#         plt.legend()
# plt.tight_layout()