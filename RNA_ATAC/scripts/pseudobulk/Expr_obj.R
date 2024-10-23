########################
##                    ##
##  Pseudotime plots  ##
##                    ##
########################

source("/data/homes/louisc/Project_Babraham/RNA/scripts/Settings.R")

args <- list()
args$sce <- "/data/homes/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
args$sce_pseudobulk <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_pseudobulk.rds"
args$sce_metacell <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/SingleCellExperiment_metacells_nodiff.rds"
args$metadata <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/clustering/PeakMatrix/sample_metadata_nodiff_after_clustering.txt.gz"
args$group_by <- "cluster"


###############
## Load data ##
###############

sce_pseudobulk <- readRDS(args$sce_pseudobulk)
sce_metacell <- readRDS(args$sce_metacell)

#####################
## Day of sampling ##
#####################

# Load cell metadata
sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  dplyr::rename(cell=sample)

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Remove cluster 6 cells
sce_nc6 <- sce[,sce$cluster!=6]

# Aggregate data by day of sampling
sce_dos <- pseudobulk_sce_fn(
  x = sce_nc6,
  assay = "counts",
  by = "batch",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_dos) <- "counts"

# Normalisation 
logcounts(sce_dos) <- log2(1e6*(sweep(counts(sce_dos),2,colSums(counts(sce_dos)),"/"))+1)

# Save
saveRDS(sce_dos, "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/SingleCellExperiment_dayofsampling.rds")

#########################
## Expression matrices ##
#####################.

if (!dir.exists("/data/homes/louisc/Project_Babraham/expr_mat")){dir.create("/data/homes/louisc/Project_Babraham/expr_mat")}

dat_mat_raw_pb <- counts(sce_pseudobulk)
dat_mat_raw_dos <- counts(sce_dos)
dat_mat_raw_mc <- counts(sce_metacell)

write.table(dat_mat_raw_pb,file="/data/homes/louisc/Project_Babraham/expr_mat/expr_raw_summed_pb.txt",col.names=T, row.names=T, sep="\t", quote=F)
write.table(dat_mat_raw_dos,file="/data/homes/louisc/Project_Babraham/expr_mat/expr_raw_summed_dos.txt",col.names=T, row.names=T, sep="\t", quote=F)
write.table(dat_mat_raw_mc,file="/data/homes/louisc/Project_Babraham/expr_mat/expr_raw_summed_mc.txt",col.names=T, row.names=T, sep="\t", quote=F)

dat_mat_norm_pb <- logcounts(sce_pseudobulk)
dat_mat_norm_dos <- logcounts(sce_dos)
dat_mat_norm_mc <- logcounts(sce_metacell)

write.table(dat_mat_norm_pb,file="/data/homes/louisc/Project_Babraham/expr_mat/expr_norm_pb.txt",col.names=T, row.names=T, sep="\t", quote=F)
write.table(dat_mat_norm_dos,file="/data/homes/louisc/Project_Babraham/expr_mat/expr_norm_dos.txt",col.names=T, row.names=T, sep="\t", quote=F)
write.table(dat_mat_norm_mc,file="/data/homes/louisc/Project_Babraham/expr_mat/expr_norm_mc.txt",col.names=T, row.names=T, sep="\t", quote=F)