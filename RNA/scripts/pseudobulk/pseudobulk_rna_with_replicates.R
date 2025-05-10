########################################
##                                    ##
##  Pseudobulk_rna_with_replicates.R  ##
##                                    ##
########################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",    help='SingleCellExperiment object')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--cell2replicate',    type="character",    help='')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--nrep',       type="integer",       default=5,      help='Number of replicates per group (cells sampled with replacement)')
p$add_argument('--min_cells',       type="integer",       default=25,      help='Minimum number of cells per replicate')
p$add_argument('--n_neighbors_clusters',     type="integer",    default=15,     help='Louvain clustering number of neighbours')
p$add_argument('--resolution',        type="double",     default=0.25,     help='Louvain clustering resolution')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--fraction_cells_per_replicate',       type="double",       default=0.3,      help='Percentage of cells per replicate')
p$add_argument('--outdir',      type="character",    help='Output directory')
p$add_argument("--incl_samples",  type="character",    help='Which samples should be included')
args <- p$parse_args(commandArgs(TRUE))

dir.create(args$outdir, showWarnings = F, recursive = T)

args$filter_differentiated <- FALSE

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  setnames(args$group_by,"group")

print(table(cell_metadata.dt$group))

# Remove groups with not enough cells
groups.to.remove <- cell_metadata.dt[,.N,by="group"] %>% .[N<=args$min_cells,group]
print(sprintf("Removing the following groups because they have less than %d cells:",args$min_cells))
print(cell_metadata.dt[group%in%groups.to.remove,.N,by="group"])
cell_metadata.dt <- cell_metadata.dt[!group%in%groups.to.remove]

##################################
## Create pseudobulk replicates ##
##################################

set.seed(args$seed1)
cell2replicate.dt <- unique(cell_metadata.dt$group) %>% map(function(i) {
  tmp <- cell_metadata.dt[group==i]
  if ((args$fraction_cells_per_replicate*nrow(tmp))<=args$min_cells) {
    ncells_per_replicate <- args$min_cells
  } else {
    ncells_per_replicate <- round(args$fraction_cells_per_replicate*nrow(tmp))
  }
  seq(1,args$nrep) %>% map(function(j) {
    tmp[sample.int(nrow(tmp),ncells_per_replicate)] %>% 
      .[,replicate:=sprintf("%s_rep%s",i,j)] %>%
      .[,c("cell","group","replicate")] %>% 
      return
  }) %>% rbindlist %>% return
}) %>% rbindlist

stats.dt <- cell2replicate.dt[,.(ncells=.N),c("group","replicate")]
print(stats.dt)
print(stats.dt[!duplicated(stats.dt$group),-2])

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(args$sce, cells=cell2replicate.dt$cell)
sce$group_by <- cell2replicate.dt$replicate

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "group_by",
  fun = "sum",
  scale = FALSE
)

assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)

##########
## Save ##
##########

# Save cell2replicate assignment
fwrite(cell2replicate.dt, file.path(args$outdir,"cell2replicate.txt.gz"), sep="\t", quote = F)

# Save SingleCellExperiment
if (args$group_by=="cluster") {
  sce_pseudobulk$celltype <- colnames(sce_pseudobulk) %>% strsplit("_rep") %>% map_chr(1)
}
saveRDS(sce_pseudobulk, file.path(args$outdir,"SingleCellExperiment_pseudobulk_with_replicates.rds"))