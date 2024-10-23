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
p.add_argument( '--incl_samples',            type=str,   required=True,       default="all",             help='Samples')
p.add_argument( '--cellranger_output',               type=str,                required=True,           help='Output from cellranger')
p.add_argument( '--outfile',               type=str,                required=True,           help='Output file (anndata)')
args = p.parse_args()

## START TEST ##
# args = {}
# # args["samples"] = ["E7.5_rep1", "E7.5_rep2", "E7.75_rep1", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2", "E8.75_rep1", "E8.75_rep2", "E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT"]
# args["samples"] = "all"
# # args["loom_files"] = [io["basedir"]+"/original/"%s/velocyto"%i for i in args["samples"]]
# args["anndata"] = "/data/homes/louisc/Project_Babraham/RNA/anndata.h5ad"
# args["metadata"] = "/data/homes/louisc/Project_Babraham/RNA/mapping/sample_metadata_after_doublets.txt.gz"
# args["cellranger_output"] = "/data/homes/louisc/Project_Babraham/results_sierra"
# args["outfile"] = "/data/homes/louisc/Project_Babraham/RNA/velocyto/anndata_velocyto.h5ad"
## END TEST ##

# convert args to dictionary
args = vars(args)

#####################
## Parse arguments ##
#####################

if not os.path.isdir(os.path.dirname(args["outfile"])):
    os.makedirs(os.path.dirname(args["outfile"]))

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
    normalise = False,
    cells = metadata.index.values
)
adata

##########################
## Load velocyto output ##
##########################

print("Loading velocyto output files...")

looms = [None for i in range(len(args["samples"]))]
for i in range(len(args["samples"])):
    print(args["samples"][i])
    args["loom_velocyto"] = "%s/%s/velocyto/%s.loom" % (args["cellranger_output"],args["samples"][i],args["samples"][i])
    looms[i] = sc.read_loom(args["loom_velocyto"], sparse=True, X_name='spliced', obs_names='CellID', obsm_names=None, var_names='Gene')
    looms[i].var_names_make_unique()
    looms[i].obs.index = looms[i].obs.index.str.replace(":","#").str.replace("x","-1")
    print(looms[i].shape)

###########################
## Parse velocyto output ##
###########################

print("Parsing velocyto output files...")

adata_loom = anndata.AnnData.concatenate(*looms, join='inner', batch_key=None, index_unique=None)
adata_loom
del looms

# Remove non-used layers to save memory
del adata_loom.layers["ambiguous"]
del adata_loom.layers["matrix"]

# Convert to sparse matrices
# from scipy.sparse import csr_matrix
# adata.X = csr_matrix(adata.X)

###########################
## Merge anndata objects ##
###########################

print("Merging original anndata with the velocyto anndata...")

adata_final = scv.utils.merge(adata, adata_loom)
adata_final

adata_final.obs.index.name = None

# TO-DO: TRANSFE .UNS AND .OBSM

# del adata_loom
# del adata

##########
## Save ##
##########

print("Saving output...")

# Logical columns need to be stored as str to avoid hdf5 complaining
# adata_final.obs["pass_rnaQC"] = adata_final.obs["pass_rnaQC"].astype("category")
adata_final.obs["pass_rnaQC"] = adata_final.obs["pass_rnaQC"].astype(str)
adata_final.obs["doublet_call"] = adata_final.obs["doublet_call"].astype(str)
# adata_final.obs["genotype"] = adata_final.obs["genotype"].astype(str)
# adata_final.obs["pass_atacQC"] = adata_final.obs["pass_atacQC"].astype(str)
# adata_final.obs = adata_final.obs.drop(["pass_rnaQC","pass_atacQC","barcode","doublet_call"], axis=1)

adata_final.write_h5ad(args["outfile"], compression="gzip")