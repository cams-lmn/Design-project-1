##########################################
##                                      ##
##  ArchR_pseudobulk_with_replicates.R  ##
##                                      ##
##########################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atac_matrix_name',    type="character",    help='ATAC matrix')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--nrep',       type="integer",       default=5,      help='Number of replicates per group (cells sampled with replacement)')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per replicate')
p$add_argument('--fraction_cells_per_replicate',       type="double",       default=0.3,      help='Percentage of cells per replicate')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

args$atac_matrix_file <- file.path(sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/archR/Matrices/%s_summarized_experiment.rds",
                                           args$atac_matrix_name))
args$outdir <- file.path(sprintf("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/archR/pseudobulk/%s/%s",
                                 args$group_by,args$atac_matrix_name))

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  setnames(args$group_by,"group")

print(cell_metadata.dt[,.N,by="group"])

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

######################
## Load ATAC matrix ##
######################

print("Loading ATAC matrix...")

atac.se <- readRDS(args$atac_matrix_file)[,cell2replicate.dt$cell] %>% as(.,"SingleCellExperiment")
atac.se$group_by <- cell2replicate.dt$replicate
assayNames(atac.se) <- "counts"

################
## Pseudobulk ##
################

print("Pseudobulking...")

atac_pseudobulk.se <- pseudobulk_sce_fn(
  x = atac.se,
  assay = "counts",
  by = "group_by",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(atac_pseudobulk.se) <- "counts"

###################
## Normalisation ##
###################

logcounts(atac_pseudobulk.se) <- log(1e6*(sweep(counts(atac_pseudobulk.se),2,colSums(counts(atac_pseudobulk.se)),"/"))+1)

# Edit colData
if (args$group_by=="celltype_genotype") {
  atac_pseudobulk.se$celltype <- colnames(atac_pseudobulk.se) %>% strsplit("-") %>% map_chr(1)
  atac_pseudobulk.se$genotype <- colnames(atac_pseudobulk.se) %>% strsplit("-") %>% map_chr(2) %>% strsplit("_rep") %>% map_chr(1)
  atac_pseudobulk.se$celltype_genotype <- sprintf("%s_%s",atac_pseudobulk.se$celltype,atac_pseudobulk.se$genotype)
} else if (args$group_by=="celltype") {
  atac_pseudobulk.se$celltype <- colnames(atac_pseudobulk.se) %>% strsplit("_rep") %>% map_chr(1)
} else if (args$group_by=="Clusters") {
  atac_pseudobulk.se$celltype <- colnames(atac_pseudobulk.se) %>% strsplit("_rep") %>% map_chr(1)
}


# Save Summarizedexperiment
saveRDS(atac_pseudobulk.se, file.path(args$outdir,"pseudobulk_with_replicates.rds"))

# Save cell2replicate assignment
fwrite(cell2replicate.dt, file.path(args$outdir,"cell2replicate.txt.gz"), sep="\t", quote = F)