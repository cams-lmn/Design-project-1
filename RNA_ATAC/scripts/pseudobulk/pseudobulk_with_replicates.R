###############################
##                           ##
##  Pseudobulk_replicates.R  ##
##                           ##
###############################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",    help='SingleCellExperiment object')
p$add_argument('--atac_matrix_names',         type="character", nargs="+",   help='ATAC matrices')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--nrep',       type="integer",       default=5,      help='Number of replicates per group (cells sampled with replacement)')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per replicate')
p$add_argument('--fraction_cells_per_replicate',       type="double",       default=0.3,      help='Percentage of cells per replicate')
p$add_argument('--outdir',      type="character",    help='Output directory')
p$add_argument('--seed1',       type="integer",       default=42,      help='Random seed1')
args <- p$parse_args(commandArgs(TRUE))

# START TEST ##
# args <- list()
# args$filter_differentiated <- TRUE
# if (args$filter_differentiated){
#   args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/sample_metadata_nodiff_after_clustering.txt.gz"
# } else {
#   args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/mofa/fast/sample_metadata_after_clustering.txt.gz"
#   
# }
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$atac_matrix_names <- c("GeneScoreMatrix_TSS","PeakMatrix","GeneScoreMatrix_distal")
# args$group_by <- "cluster"
# args$nrep <- 5
# args$min_cells <- 25
# args$fraction_cells_per_replicate <- 0.30
# args$seed1 <- 42
# args$outdir <- file.path(sprintf("/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/%s",args$group_by))
# END TEST ##

#####################
## Define settings ##
#####################

dir.create(args$outdir, showWarnings = F, recursive = T)

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(args$metadata) %>%
  # .[,celltype_genotype:=as.character(NA)] %>% .[!is.na(genotype),celltype_genotype:=sprintf("%s-%s",celltype,genotype)] %>%
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
      .[,c("sample","group","replicate")] %>% 
      return
  }) %>% rbindlist %>% return
}) %>% rbindlist

stats.dt <- cell2replicate.dt[,.(ncells=.N),c("group","replicate")]
print(stats.dt)
print(stats.dt[!duplicated(stats.dt$group),-2])

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(args$sce, cells=cell2replicate.dt$sample)
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
sce_pseudobulk$celltype <- colnames(sce_pseudobulk) %>% strsplit("_rep") %>% map_chr(1)
dir.create(paste0(args$outdir,"/RNA"), showWarnings = F, recursive = T)
saveRDS(sce_pseudobulk, file.path(args$outdir,"RNA/pseudobulk_with_replicates.rds"))


for (i in 1:length(args$atac_matrix_names)){
  args$atac_matrix_name <- args$atac_matrix_names[i]
  args$atac_matrix_file <- file.path(sprintf("/data/homes/louisc/Project_Babraham/ATAC/archR/Matrices/%s_summarized_experiment.rds",
                                             args$atac_matrix_name))
  
  ######################
  ## Load ATAC matrix ##
  ######################
  
  atac.se <- readRDS(args$atac_matrix_file)[,cell2replicate.dt$sample] %>% as(.,"SingleCellExperiment")
  atac.se$group_by <- cell2replicate.dt$replicate
  assayNames(atac.se) <- "counts"
  
  ################
  ## Pseudobulk ##
  ################
  
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
  
  ##########
  ## Save ##
  ##########
  
  # Save Summarizedexperiment
  atac_pseudobulk.se$celltype <- colnames(atac_pseudobulk.se) %>% strsplit("_rep") %>% map_chr(1)
  dir.create(paste0(args$outdir,"/",args$atac_matrix_name), showWarnings = F, recursive = T)
  saveRDS(atac_pseudobulk.se, sprintf("%s/%s/pseudobulk_with_replicates.rds",args$outdir,args$atac_matrix_name))
}

