###################
##               ##
##  CellPhoneDb  ##
##               ##
###################

######################
## Import libraries ##
######################

import os
# from re import search
# from dfply import *
import sys
import matplotlib.pyplot as plt
import ktplotspy as kpy
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as s
import re
from copy import copy
# import math
# import seaborn as sns
# import importlib
import argparse
from cellphonedb.utils import db_utils, db_releases_utils
from IPython.display import HTML, display
from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
import anndata

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

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--cpdb_version',               type=str,                required=True,           help='CellPhoneDb version')
p.add_argument( '--db_path',               type=str,                required=True,           help='CellPhoneDb path')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file (input)')
p.add_argument( '--metadata_out',               type=str,                required=True,           help='Cell metadata file (output)')
p.add_argument( '--incl_samples',            type=str,    required=True,       default="all",             help='Samples')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file (input)')
p.add_argument( '--anndata_out',               type=str,                required=True,           help='Anndata file (output)')
p.add_argument( '--degs',               type=str,                required=True,           help='Differentially expressed genes file (input)')
p.add_argument( '--degs_out',               type=str,                required=True,           help='Differentially expressed genes file (output)')
p.add_argument( '--celltypes_to_use',               type=int,                required=True, nargs="+",            help='Celltypes to use')
p.add_argument( '--gene_input',               type=str,                required=True,           help='Genes (input)')
p.add_argument( '--complex_input',               type=str,                required=True,           help='Gene complexes (input)')
args = p.parse_args()

# ## START TEST ##
# args = {}
# args["cpdb_version"] = "v5.0.0"
# args["db_path"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/db_input_files"
# args["metadata"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
# args["metadata_out"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/metadata.tsv"
# args["incl_samples"] = "nodiff"
# args["outdir"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/"
# args["anndata"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/anndata_nodiff.h5ad"
# args["anndata_out"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/normalised_log_counts.h5ad"
# args["degs"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
# args["degs_out"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/degs.tsv"
# args["celltypes_to_use"] = [ 1, 2, 3, 4, 5 ]
# args["gene_input"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/db_input_files/v5.0.0/gene_input.csv"
# args["complex_input"] = "/data/homes/louisc/Project_Babraham/RNA_ATAC/CellPhoneDb/db_input_files/v5.0.0/complex_input.csv"
# ## END TEST ##

# convert args to dictionary
args = vars(args)

opts = {}
opts["celltype_colors"] = ["#DAF7A6","#FFC300","#DDAA05","#C70039","#900C3F","#581845"]

#####################
## Parse arguments ##
#####################

print(args)

if not os.path.isdir(args["outdir"]):
    args["outdir"] = Path(args["outdir"])

if args["incl_samples"]=="all":
    args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L"]
elif args["incl_samples"]=="nodiff":
    args["samples"] = ["d0","d1","d3","d5","d7","d10","d14","d18"]

#######################
## Download database ##
#######################

# Version of the databse
print(args["cpdb_version"])

# Path where the input files to generate the database are located
cpdb_target_dir = os.path.join(args["db_path"], args["cpdb_version"])
db_utils.download_database(cpdb_target_dir, args["cpdb_version"])
args["cpdb_file_path"] = args["db_path"] + "/" + args["cpdb_version"] + "/cellphonedb.zip"

###################
## Load metadata ##
###################

metadata = pd.read_table(args["metadata"])[["sample","cluster"]].rename(columns={"sample":"barcode_sample","cluster":"cell_type"})
metadata.head(3)
metadata.shape
metadata.to_csv(args["metadata_out"], index=False, sep="\t") 
#metadata =  pd.read_csv(args["metadata_out"], sep='\t')

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

sum(np.array(list(adata.obs.index)) == np.array(list(metadata['barcode_sample'])))

adata.obs.index.names = ["barcode_sample"]
adata.obs = adata.obs[["nFeature_RNA","nCount_RNA"]]
adata.obs = adata.obs.rename(columns={"nFeature_RNA":"n_genes","nCount_RNA":"n_counts"}) 
adata.obs["cell_labels"] = list(metadata['cell_type'])
adata.var = adata.var.rename(columns={"gene":"gene_ids"})
adata.var["feature_types"] = "Gene Expression"

del adata.layers
del adata.uns

adata.write_h5ad(args["anndata_out"])

#adata = anndata.read_h5ad(args["anndata_out"])

#################################
## Load and transform DEG info ##
#################################.

deg_data = pd.read_table(args["degs"])
for i in range(1,6):
    print(i)
    # Marker genes cluster i
    # deg_i_index = deg_data.index[deg_data["clust" + str(i) + "_spec"]==True]
    # if i == 1:
    #     deg_data_transf = [pd.concat([pd.DataFrame([str(i)]*len(deg_i_index),index=deg_i_index), deg_data["Gene"][deg_data.index.isin(deg_i_index)]], axis=1, ignore_index=True)]
    # else :
    #     deg_data_transf.append(pd.concat([pd.DataFrame([str(i)]*len(deg_i_index),index=deg_i_index), deg_data["Gene"][deg_data.index.isin(deg_i_index)]], axis=1, ignore_index=True))
    # DE genes cluster i
    first_iteration = True
    for j in range(1,6):
        # print(j)
        if j != i:
            if j > i:
                deg_data_sub = deg_data[~pd.isnull(deg_data["Sign" + str(i) + "vs" + str(j)])]
                if first_iteration:
                    deg_sub_i_index = deg_data_sub.index[np.array(deg_data_sub["Sign" + str(i) + "vs" + str(j)].tolist()) & np.array((deg_data_sub["logFC_" + str(i) + "vs" + str(j)]>0).tolist())].tolist()
                else:
                    deg_sub_i_index = deg_sub_i_index.append(deg_data_sub.index[np.array(deg_data_sub["Sign" + str(i) + "vs" + str(j)].tolist()) & np.array((deg_data_sub["logFC_" + str(i) + "vs" + str(j)]>0).tolist())].tolist()).unique().sort()
            else:
                deg_data_sub = deg_data[~pd.isnull(deg_data["Sign" + str(j) + "vs" + str(i)])]
                if first_iteration:
                    deg_sub_i_index = deg_data_sub.index[np.array(deg_data_sub["Sign" + str(j) + "vs" + str(i)].tolist()) & np.array((deg_data_sub["logFC_" + str(j) + "vs" + str(i)]<0).tolist())].tolist()
                else: 
                    deg_sub_i_index = deg_sub_i_index.append(deg_data_sub.index[np.array(deg_data_sub["Sign" + str(j) + "vs" + str(i)].tolist()) & np.array((deg_data_sub["logFC_" + str(j) + "vs" + str(i)]<0).tolist())].tolist()).unique().sort()
    if i == 1:
        deg_data_transf = [pd.concat([pd.DataFrame([str(i)]*len(deg_sub_i_index), index=deg_sub_i_index),deg_data_sub["Gene"][deg_sub_i_index],deg_data_sub["gene_cluster"][deg_sub_i_index]], axis=1)]
    else:
        deg_data_transf.append(pd.concat([pd.DataFrame([str(i)]*len(deg_sub_i_index), index=deg_sub_i_index),deg_data_sub["Gene"][deg_sub_i_index],deg_data_sub["gene_cluster"][deg_sub_i_index]], axis=1))

deg_data_transf = pd.concat(deg_data_transf, ignore_index=True)
deg_data_transf.columns = ["cell_type","gene","gene_cluster"]
deg_data_transf.head()
deg_data_transf.shape

deg_data_transf.to_csv(args["degs_out"], index=False, sep="\t") 

############################
## Load Microenvironments ##
############################

# microenv = pd.read_csv(microenvs_file_path,
#                        sep = '\t')
# microenv.head(3)

# microenv.groupby('microenvironment', group_keys = False)['cell_type'] \
#     .apply(lambda x : list(x.value_counts().index))


#     from cellphonedb.src.core.methods import cpdb_analysis_method

########################
## Run basic analysis ##
########################

cpdb_results_basic = cpdb_analysis_method.call(
    cpdb_file_path = args["cpdb_file_path"],   # mandatory: CellphoneDB database zip file.
    meta_file_path = args["metadata_out"],     # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = args["anndata_out"],    # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',               # defines the gene annotation in counts matrix.
    score_interactions = True,                 # optional: whether to score interactions or not. 
    output_path = args["outdir"],              # Path to save results    microenvs_file_path = None,
    separator = '|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    threads = 5,                               # number of threads to use in the analysis.
    threshold = 0.1,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                      # Sets the rounding for the mean values in significan_means.
    debug = False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_suffix = args["incl_samples"]       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
)

cn_oi = ["%s|%s" % (cluster1, cluster2) for cluster1, cluster2 in [(args["celltypes_to_use"][i], args["celltypes_to_use"][i]) for i in range(0, len(args["celltypes_to_use"]))]]
cpdb_results_basic['means_result'].head(2)[cn_oi]
cpdb_results_basic['interaction_scores'].head(2)[cn_oi]

#######################
## Run stat analysis ##
#######################

cpdb_results_stat = cpdb_statistical_analysis_method.call(
    cpdb_file_path = args["cpdb_file_path"],         # mandatory: CellphoneDB database zip file.
    meta_file_path = args["metadata_out"],           # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = args["anndata_out"],          # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 5,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = args["outdir"],                    # Path to save results    microenvs_file_path = None,
    output_suffix = args["incl_samples"]             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )

cn_oi = ["%s|%s" % (cluster1, cluster2) for cluster1, cluster2 in [(args["celltypes_to_use"][i], args["celltypes_to_use"][i]) for i in range(0, len(args["celltypes_to_use"]))]]
cpdb_results_stat['pvalues'].head(2)[cn_oi]
cpdb_results_stat['significant_means'].head(2)[["id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted",
                                                'annotation_strategy', 'is_integrin', 'directionality', 'classification', 'rank'] + cn_oi]

# Number of significant interactions
p_c_hm = kpy.plot_cpdb_heatmap(pvals = cpdb_results_stat['pvalues'],
                      degs_analysis = False,
                      figsize = (5, 5),
                      title = "Sum of significant interactions")
p_c_hm.savefig(args["outdir"]+"/heatmap_sign_interactions")

# Extract all significant interactions
cn_info = ["id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", 'annotation_strategy', 'is_integrin', 'directionality', 'classification', 'rank']
for cn in cn_oi:
    print(cn)
    cn_index = cpdb_results_stat['significant_means'].index[~np.isnan(cpdb_results_stat['significant_means'][cn]).values]
    if cn == "1|1":
        res_final_stat = [pd.concat([pd.DataFrame([cn]*sum(cpdb_results_stat["significant_means"].index.isin(cn_index)),index=cn_index),
                                     cpdb_results_stat["significant_means"][cpdb_results_stat["significant_means"].index.isin(cn_index)][cn_info + [cn]],
                                     cpdb_results_stat["interaction_scores"][cpdb_results_stat["interaction_scores"].index.isin(cn_index)][cn],
                                     cpdb_results_stat["pvalues"][cpdb_results_stat["pvalues"].index.isin(cn_index)][cn]], axis=1, ignore_index=True)]
    else :
        res_final_stat.append(pd.concat([pd.DataFrame([cn]*sum(cpdb_results_stat["significant_means"].index.isin(cn_index)),index=cn_index),
                                         cpdb_results_stat["significant_means"][cpdb_results_stat["significant_means"].index.isin(cn_index)][cn_info + [cn]],
                                         cpdb_results_stat["interaction_scores"][cpdb_results_stat["interaction_scores"].index.isin(cn_index)][cn],
                                         cpdb_results_stat["pvalues"][cpdb_results_stat["pvalues"].index.isin(cn_index)][cn]], axis=1, ignore_index=True))

res_final_stat = pd.concat(res_final_stat, ignore_index=True)
res_final_stat.columns = ["Comparison"] + cn_info + ["significant_mean","interaction_score","pvalue"]

# Make deconv & deconv_pct object for indexing later
#deconv = cpdb_results_stat["deconvoluted"]
#deconv.index = deconv["gene"]
#deconv.loc[np.nan] = [np.nan]*len(deconv.iloc[0])
# deconv = deconv[["1","2","3","4","5"]].loc[~deconv.index.duplicated()]

deconv = cpdb_results_stat["deconvoluted"]
deconv.index = deconv["id_cp_interaction"]

# deconv_pct = cpdb_results_stat["deconvoluted_percents"]
# deconv_pct.index = deconv_pct["gene"]
# deconv_pct.loc[np.nan] = [np.nan]*len(deconv_pct.iloc[0])
# deconv_pct = deconv_pct[["1","2","3","4","5"]].loc[~deconv_pct.index.duplicated()]

deconv_pct = cpdb_results_stat["deconvoluted_percents"]
deconv_pct.index = deconv_pct["id_cp_interaction"]

# # Parse complexes and make means & pct object based on min value for all partners (add them to previous deconv & conv object)
# gene_input = pd.read_csv(args["gene_input"])
# gene_input.index=gene_input["uniprot"]
# complex_input = pd.read_csv(args["complex_input"])
# complex_input.index = complex_input["complex_name"]

# for i in complex_input["complex_name"]:
#     prot_i = complex_input[['uniprot_1', 'uniprot_2', 'uniprot_3', 'uniprot_4', 'uniprot_5']].loc[i]
#     prot_i = prot_i[~pd.isnull(prot_i)].tolist()
#     gene_i = gene_input["gene_name"].loc[prot_i].tolist()
#     if sum(deconv.index.isin(gene_i))>0:
#         gene_i = pd.DataFrame(gene_i)[0][pd.DataFrame(gene_i).isin(deconv.index)[0]].tolist()
#         if len(deconv.loc[gene_i].index)>1:
#             deconv.loc[i] = deconv.loc[gene_i].min()
#             deconv_pct.loc[i] = deconv_pct.loc[gene_i].min()
#         elif len(deconv.loc[gene_i].index)==1:
#             deconv.loc[i] = deconv.loc[gene_i].iloc[0]
#             deconv_pct.loc[i] = deconv_pct.loc[gene_i].iloc[0]

# Change res_final_stat to display complex names as gene when nan
new_gene_a = [""]*len(res_final_stat["gene_a"])
new_gene_b = [""]*len(res_final_stat["gene_b"])
for i in range(0,len(new_gene_a)):
    unique_deconv_i = np.unique(deconv.loc[res_final_stat["id_cp_interaction"][i]]["1"])
    if type(deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"])==str:
        new_gene_a[i] = deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"]
        new_gene_b[i] = deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"]
    elif len(deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"])==2:
        new_gene_a[i] = deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"][0]
        new_gene_b[i] = deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"][1]
    elif len(unique_deconv_i)==2:
        new_gene_a[i] = "_".join(deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"][deconv.loc[res_final_stat["id_cp_interaction"][i]]["1"]==unique_deconv_i[0]])
        new_gene_b[i] = "_".join(deconv.loc[res_final_stat["id_cp_interaction"][i]]["gene"][deconv.loc[res_final_stat["id_cp_interaction"][i]]["1"]==unique_deconv_i[1]])
    else: 
        new_gene_a[i] = res_final_stat["gene_a"][i]
        new_gene_b[i] = res_final_stat["gene_b"][i]

res_final_stat["gene_a"] = new_gene_a
res_final_stat["gene_b"] = new_gene_b

deconv.index = deconv["id_cp_interaction"] + "_" + deconv["gene"]
deconv.loc[np.nan] = [np.nan]*len(deconv.iloc[0])
deconv = deconv[["1","2","3","4","5"]].loc[~deconv.index.duplicated()]

deconv_pct.index = deconv_pct["id_cp_interaction"] + "_" + deconv_pct["gene"]
deconv_pct.loc[np.nan] = [np.nan]*len(deconv_pct.iloc[0])
deconv_pct = deconv_pct[["1","2","3","4","5"]].loc[~deconv_pct.index.duplicated()]

#new_gene_a = res_final_stat["partner_a"][pd.isnull(res_final_stat["gene_a"])].str.replace("^complex:","")
#res_final_stat["gene_a"][pd.isnull(res_final_stat["gene_a"])] = new_gene_a
    
#new_gene_b = res_final_stat["partner_b"][pd.isnull(res_final_stat["gene_b"])].str.replace("^complex:","")
#res_final_stat["gene_b"][pd.isnull(res_final_stat["gene_b"])] = new_gene_b

# Select res_final_stat data from decon & deconv_pct objects
# deconv_gene_a = deconv.loc[res_final_stat["gene_a"]]
# deconv_gene_a.columns = [i + j for i, j in zip(["gene_a_"]*5,deconv_gene_a.columns)]
# deconv_gene_b = deconv.loc[res_final_stat["gene_b"]]
# deconv_gene_b.columns = [i + j for i, j in zip(["gene_b_"]*5,deconv_gene_b.columns)]

# deconv_pct_gene_a = deconv_pct.loc[res_final_stat["gene_a"]]
# deconv_pct_gene_a.columns = [i + j for i, j in zip(["gene_a_pct_"]*5,deconv_pct_gene_a.columns)]
# deconv_pct_gene_b = deconv_pct.loc[res_final_stat["gene_b"]]
# deconv_pct_gene_b.columns = [i + j for i, j in zip(["gene_b_pct_"]*5,deconv_pct_gene_b.columns)]

deconv_gene_a = copy(deconv.iloc[range(0,len(res_final_stat)),:])
deconv_gene_a.loc[:] = np.nan
deconv_gene_a.index = range(0,len(deconv_gene_a.index))
mat_gene_a = deconv.loc[(res_final_stat["id_cp_interaction"][pd.notnull(res_final_stat["gene_a"])] + "_" + [re.sub("_[A-Z,0-9,_]+$","",x) for x in res_final_stat["gene_a"][pd.notnull(res_final_stat["gene_a"])]])]
mat_gene_a.index = deconv_gene_a.index[pd.notnull(res_final_stat["gene_a"]).tolist()]
deconv_gene_a.loc[pd.notnull(res_final_stat["gene_a"]).tolist(),:] = mat_gene_a
deconv_gene_a.columns = [i + j for i, j in zip(["gene_a_"]*5,deconv_gene_a.columns)]

deconv_gene_b = copy(deconv.iloc[range(0,len(res_final_stat)),:])
deconv_gene_b.loc[:] = np.nan
deconv_gene_b.index = range(0,len(deconv_gene_b.index))
mat_gene_b = deconv.loc[(res_final_stat["id_cp_interaction"][pd.notnull(res_final_stat["gene_b"])] + "_" + [re.sub("_[A-Z,0-9,_]+$","",x) for x in res_final_stat["gene_b"][pd.notnull(res_final_stat["gene_b"])]])]
mat_gene_b.index = deconv_gene_a.index[pd.notnull(res_final_stat["gene_b"]).tolist()]
deconv_gene_b.loc[pd.notnull(res_final_stat["gene_b"]).tolist(),:] = mat_gene_b
deconv_gene_b.columns = [i + j for i, j in zip(["gene_b_"]*5,deconv_gene_b.columns)]

deconv_pct_gene_a = copy(deconv_pct.iloc[range(0,len(res_final_stat)),:])
deconv_pct_gene_a.loc[:] = np.nan
deconv_pct_gene_a.index = range(0,len(deconv_pct_gene_a.index))
mat_pct_gene_a = deconv_pct.loc[(res_final_stat["id_cp_interaction"][pd.notnull(res_final_stat["gene_a"])] + "_" + [re.sub("_[A-Z,0-9,_]+$","",x) for x in res_final_stat["gene_a"][pd.notnull(res_final_stat["gene_a"])]])]
mat_pct_gene_a.index = deconv_pct_gene_a.index[pd.notnull(res_final_stat["gene_a"]).tolist()]
deconv_pct_gene_a.loc[pd.notnull(res_final_stat["gene_a"]).tolist(),:] = mat_pct_gene_a
deconv_pct_gene_a.columns = [i + j for i, j in zip(["gene_a_pct_"]*5,deconv_pct_gene_a.columns)]

deconv_pct_gene_b = copy(deconv_pct.iloc[range(0,len(res_final_stat)),:])
deconv_pct_gene_b.loc[:] = np.nan
deconv_pct_gene_b.index = range(0,len(deconv_pct_gene_b.index))
mat_pct_gene_b = deconv_pct.loc[(res_final_stat["id_cp_interaction"][pd.notnull(res_final_stat["gene_b"])] + "_" + [re.sub("_[A-Z,0-9,_]+$","",x) for x in res_final_stat["gene_b"][pd.notnull(res_final_stat["gene_b"])]])]
mat_pct_gene_b.index = deconv_pct_gene_b.index[pd.notnull(res_final_stat["gene_b"]).tolist()]
deconv_pct_gene_b.loc[pd.notnull(res_final_stat["gene_b"]).tolist(),:] = mat_pct_gene_b
deconv_pct_gene_b.columns = [i + j for i, j in zip(["gene_b_pct_"]*5,deconv_pct_gene_b.columns)]

for i in range(0,len(res_final_stat.index)):
    if pd.isnull(res_final_stat["gene_a"][i]):
        list_to_add = [np.nan]
        if pd.isnull(res_final_stat["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(','.join(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist()))
    else: 
        if len(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_stat["gene_a"][i]].tolist())==0:
            list_to_add = [np.nan]
        else: 
            list_to_add = [','.join(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_stat["gene_a"][i]].tolist())]
        if pd.isnull(res_final_stat["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(','.join(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist()))
    if i == 0:
        deg_stat = [list_to_add]
    else: 
        deg_stat.append(list_to_add)       

deg_stat = pd.DataFrame(deg_stat)
deg_stat.columns = ["deg_gene_a","deg_gene_b"]    

for i in range(0,len(res_final_stat.index)):
    if pd.isnull(res_final_stat["gene_a"][i]):
        list_to_add = [np.nan]
        if pd.isnull(res_final_stat["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist()[0])
    else: 
        if len(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_stat["gene_a"][i]].tolist())==0:
            list_to_add = [np.nan]
        else: 
            list_to_add = [deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_stat["gene_a"][i]].tolist()[0]]
        if pd.isnull(res_final_stat["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_stat["gene_b"][i]].tolist()[0])
    if i == 0:
        gc_stat = [list_to_add]
    else: 
        gc_stat.append(list_to_add)       

gc_stat = pd.DataFrame(gc_stat)
gc_stat.columns = ["gc_gene_a","gc_gene_b"] 

res_final_stat_elab = pd.concat([res_final_stat, deconv_gene_a, deconv_gene_b, deconv_pct_gene_a, deconv_pct_gene_b, deg_stat, gc_stat],axis=1)
res_final_stat_elab.to_csv(args["outdir"] + "/CellPhone_res_stat_" + args["cpdb_version"] + "_" + args["incl_samples"] + ".txt", sep="\t", index=False)

######################
## Run deg analysis ##
######################

cpdb_results_deg = cpdb_degs_analysis_method.call(
    cpdb_file_path = args["cpdb_file_path"],                    # mandatory: CellphoneDB database zip file.
    meta_file_path = args["metadata_out"],                      # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = args["anndata_out"],                     # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    degs_file_path = args["degs_out"],                          # mandatory: tsv file with DEG to account.
    counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
                                                                # optional: defines cell types and their active TFs.
                                                                # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                                  # optional: whether to score interactions or not. 
    threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
    separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_path = args["outdir"],                               # Path to save results
    output_suffix = args["incl_samples"],                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
    threads = 25
    )

# Number of relevant interactions
p_c_hm = kpy.plot_cpdb_heatmap(pvals = cpdb_results_deg['relevant_interactions'],
                      degs_analysis=True,
                      figsize=(5, 5),
                      title="Sum of significant interactions")
p_c_hm.savefig(args["outdir"]+"/heatmap_relevant_interactions")

# Extract all significant interactions
cn_info = ["id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", 'annotation_strategy', 'is_integrin', 'directionality', 'classification', 'rank']
for cn in cn_oi:
    print(cn)
    cn_index = cpdb_results_deg['relevant_interactions'].index[cpdb_results_deg['relevant_interactions'][cn]==1]
    if cn == "1|1":
        res_final_deg = [pd.concat([pd.DataFrame([cn]*sum(cpdb_results_deg["significant_means"].index.isin(cn_index)),index=cn_index),
                                    cpdb_results_deg["significant_means"][cpdb_results_deg["significant_means"].index.isin(cn_index)][cn_info + [cn]],
                                    cpdb_results_deg["interaction_scores"][cpdb_results_deg["interaction_scores"].index.isin(cn_index)][cn]], axis=1, ignore_index=True)]
    else :
        res_final_deg.append(pd.concat([pd.DataFrame([cn]*sum(cpdb_results_deg["significant_means"].index.isin(cn_index)),index=cn_index),
                                        cpdb_results_deg["significant_means"][cpdb_results_deg["significant_means"].index.isin(cn_index)][cn_info + [cn]],
                                        cpdb_results_deg["interaction_scores"][cpdb_results_deg["interaction_scores"].index.isin(cn_index)][cn]], axis=1, ignore_index=True))

res_final_deg = pd.concat(res_final_deg, ignore_index=True)
res_final_deg.columns = ["Comparison"] + cn_info + ["significant_mean","interaction_score"]

for i in range(0,len(res_final_deg.index)):
    #print(i)
    if pd.isnull(res_final_deg["gene_a"][i]):
        if pd.isnull(res_final_deg["gene_b"][i]):
            if i == 0:
                pct_expr_deg = [[np.nan,np.nan,np.nan,np.nan]]
            else:
                pct_expr_deg.append([np.nan,np.nan,np.nan,np.nan])
        else:
            if i == 0:
                pct_expr_deg = [[np.nan,
                                  cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                  np.nan,
                                  cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0]]]
            else:
                pct_expr_deg.append([np.nan,
                                      cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                      np.nan,
                                      cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0]])    
    else: 
        if pd.isnull(res_final_deg["gene_b"][i]):
            if i == 0:
                pct_expr_deg = [[cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                  np.nan,
                                  cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                  np.nan]]
            else:
                pct_expr_deg.append([cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                      np.nan,
                                      cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                      np.nan])
        else: 
            if i == 0:
                pct_expr_deg = [[cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                  cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                  cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                  cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0]]]
            else:
                pct_expr_deg.append([cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                      cpdb_results_deg["deconvoluted_percents"][cpdb_results_deg["deconvoluted_percents"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                      cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_a"][i]][res_final_deg["Comparison"][i][0]].tolist()[0],
                                      cpdb_results_deg["deconvoluted"][cpdb_results_deg["deconvoluted"]["gene_name"]==res_final_deg["gene_b"][i]][res_final_deg["Comparison"][i][0]].tolist()[0]])

pct_expr_deg = pd.DataFrame(pct_expr_deg)
pct_expr_deg.columns = ["pct_gene_a","pct_gene_b","exprs_gene_a","exprs_gene_b"]

for i in range(0,len(res_final_deg.index)):
    if pd.isnull(res_final_deg["gene_a"][i]):
        list_to_add = [np.nan]
        if pd.isnull(res_final_deg["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(','.join(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist()))
    else: 
        if len(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_deg["gene_a"][i]].tolist())==0:
            list_to_add = [np.nan]
        else: 
            list_to_add = [','.join(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_deg["gene_a"][i]].tolist())]
        if pd.isnull(res_final_deg["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(','.join(deg_data_transf["cell_type"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist()))
    if i == 0:
        deg_deg = [list_to_add]
    else: 
        deg_deg.append(list_to_add)       

deg_deg = pd.DataFrame(deg_deg)
deg_deg.columns = ["deg_gene_a","deg_gene_b"]    

for i in range(0,len(res_final_deg.index)):
    if pd.isnull(res_final_deg["gene_a"][i]):
        list_to_add = [np.nan]
        if pd.isnull(res_final_deg["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist()[0])
    else: 
        if len(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_deg["gene_a"][i]].tolist())==0:
            list_to_add = [np.nan]
        else: 
            list_to_add = [deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_deg["gene_a"][i]].tolist()[0]]
        if pd.isnull(res_final_deg["gene_b"][i]):
            list_to_add.append(np.nan)
        else:
            if len(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist())==0:
                list_to_add.append(np.nan)
            else: 
                list_to_add.append(deg_data_transf["gene_cluster"][deg_data_transf["gene"]==res_final_deg["gene_b"][i]].tolist()[0])
    if i == 0:
        gc_deg = [list_to_add]
    else: 
        gc_deg.append(list_to_add)       

gc_deg = pd.DataFrame(gc_deg)
gc_deg.columns = ["gc_gene_a","gc_gene_b"]

res_final_deg_elab = pd.concat([res_final_deg, pct_expr_deg,deg_deg,gc_deg],axis=1)
res_final_deg_elab.to_csv(args["outdir"] + "/CellPhone_res_deg_" + args["cpdb_version"] + "_" + args["incl_samples"] + ".txt", sep="\t", index=False)

# Comparison res_final_deg & res_final_stat
res_final_deg[~res_final_deg["id_cp_interaction"].isin(res_final_stat["id_cp_interaction"]).values]
res_final_stat[~res_final_stat["id_cp_interaction"].isin(res_final_deg["id_cp_interaction"]).values]