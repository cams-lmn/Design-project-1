# Hwanseok Jeong
# From what Fiebe already did, I modified so that the code can be run for all the samples.
# Also changed so that depthcor and UMAP plots can be saved.

library(dplyr)
library(data.table)
library(Seurat)
library(GenomicRanges)
library(Signac)
library(ggplot2)
library(ggpubr)


#####################
## Define settings ##
#####################

# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts <- list()
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "#"

##############################
## Load and merge data sets ##
##############################


# List of directories and corresponding names
directories <- c("C1", "C2", "C3", "Early1", "Early2", "Early3")
#base_path <- "/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/processed_GEO_data"
base_path <- "C:/Users/Hwanseok Jeong/Desktop/1-1/Design project/processed_GEO_data"

# Loop through each directory
for (i in directories) {
  # Construct the full path
  dir_path <- file.path(base_path, i)
  
  # Set the working directory
  setwd(dir_path)
  cat("Working directory set to:", getwd(), "\n")
  
  count_mtx <- list()
  cell.info <- list()
  ATAC.info <- list()
  
  # Load ATAC metadata
  ATAC.loc <- sprintf("atac_features.tsv.gz")
  ATAC.info[[i]] <- fread(ATAC.loc, header = FALSE, select = c(1, 4, 5, 6)) %>%
    setnames(c("peak", "chromosome", "start", "stop")) %>%
    .[, idx := 1:.N]
  dim(ATAC.info[[i]])
  
  # Load cell metadata
  barcode.loc <- sprintf("barcodes.tsv.gz")
  cell.info[[i]] <- fread(barcode.loc, header=F) %>%
    setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  
  # Load matrix
  matrix.loc <- sprintf("atac_matrix.mtx.gz")
  count_mtx[[i]] <- Matrix::readMM(matrix.loc)[ATAC.info[[i]]$idx,]
  stopifnot(nrow(cell.info[[i]])==ncol(count_mtx[[i]]))
  rownames(count_mtx[[i]]) <- ATAC.info[[i]]$peak
  colnames(count_mtx[[i]]) <- cell.info[[i]]$cell
  
  # Basic filtering
  cell.info[[i]] <- cell.info[[i]][Matrix::colSums(count_mtx[[i]]) <= 1.5e4, ]
  count_mtx[[i]] <- count_mtx[[i]][, Matrix::colSums(count_mtx[[i]]) <= 1.5e4]
  
  dim_count_mtx <- lapply(count_mtx,dim)
  dim_count_mtx <- do.call(cbind,dim_count_mtx)
  rownames(dim_count_mtx) <- c("# features","# cells")
  print(dim_count_mtx)
  
  #######################
  ## Keep common peaks ##
  #######################
  
  ATAC_peaks <- Reduce("intersect",lapply(count_mtx,rownames))
  for (j in 1:length(count_mtx)) {
    count_mtx[[j]] <- count_mtx[[j]][ATAC_peaks,]
  }
  
  stopifnot(length(unique(lapply(count_mtx,nrow)))==1)
  stopifnot(length(unique(lapply(count_mtx,rownames)))==1)
  
  #################
  ## Concatenate ##
  #################
  
  # Concatenate cell metadata
  cell.info <- rbindlist(cell.info)
  rownames(cell.info) <- cell.info$cell
  
  # Concatenate matrices
  count_mtx <- do.call("cbind",count_mtx)
  colnames(count_mtx) <- cell.info$cell
  
  ##################
  ## Filter peaks ##
  ##################
  
  # Remove duplicated peaks
  count_mtx <- count_mtx[!duplicated(rownames(count_mtx)),]
  
  # Sanity checks
  stopifnot(sum(duplicated(rownames(count_mtx)))==0)
  stopifnot(sum(duplicated(colnames(count_mtx)))==0)
  
  ##########################
  ## Create Seurat object ##
  ##########################
  
  cell.info.to.seurat <- cell.info[cell%in%colnames(count_mtx)] %>% setkey(cell) %>% .[colnames(count_mtx)] %>% as.data.frame
  rownames(cell.info.to.seurat) <- cell.info.to.seurat$cell
  stopifnot(rownames(cell.info.to.seurat)==colnames(count_mtx))
  stopifnot(sum(is.na(rownames(cell.info.to.seurat$cell)))==0)
  
  # Create a GRanges object from the data
  granges <- GRanges(
    seqnames = ATAC.info[[i]]$chromosome,
    ranges = IRanges(start = ATAC.info[[i]]$start, end = ATAC.info[[i]]$stop)
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = count_mtx,
    ranges = granges,     
    min.features = 200
  )
  
  seurat <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = cell.info.to.seurat
  )
  
  head(seurat@meta.data)
  
  metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
    .[,c("cell","barcode","sample","nFeature_ATAC","nCount_ATAC")]
  
  fwrite(metadata, "metadata.txt.gz", quote=F, na="NA", sep="\t")
  # fwrite(cell.info, paste0(args$outdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
  # fwrite(gene.info, paste0(args$outdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
  
  ## Save seurat file ##
  
  #setwd("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC")
  setwd("C:/Users/Hwanseok Jeong/Desktop/1-1/Design project")
  if (!file.exists("seurat")) {
    dir.create("seurat")
  }
  #setwd("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/seurat")
  setwd("C:/Users/Hwanseok Jeong/Desktop/1-1/Design project/seurat")

  filename <- paste0("seurat_", i, ".rds")
  saveRDS(seurat, filename)
  cat("Seurat object saved as:", filename, "\n")
  
  
  
  pbmc <- RunTFIDF(seurat)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  pbmc <- RunSVD(pbmc)
  
  
  
  ## Save DepthCor plot
  
  #setwd("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC")
  setwd("C:/Users/Hwanseok Jeong/Desktop/1-1/Design project")
  if (!file.exists("depthcor")) {
    dir.create("depthcor")
  }
  setwd("C:/Users/Hwanseok Jeong/Desktop/1-1/Design project/depthcor")
  #setwd("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/depthcor")
  
  depthcor_plot_filename <- paste0("depthcor_", i, ".pdf")
  pdf(depthcor_plot_filename)
  print(DepthCor(pbmc))
  dev.off()
  cat("DepthCor plot saved as:", depthcor_plot_filename, "\n")
  
  
  

  pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = c(2:4, 1:10))
  pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 1:15)
  pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
  DimPlot(object = pbmc, label = TRUE) + NoLegend()
  
  
  
  ## Save UMAP plot
  
  #setwd("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC")
  setwd("C:/Users/Hwanseok Jeong/Desktop/1-1/Design project")
  if (!file.exists("UMAP")) {
    dir.create("UMAP")
  }
  setwd("C:/Users/Hwanseok Jeong/Desktop/1-1/Design project/UMAP")
  #setwd("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/UMAP")
  
  umap_plot <- DimPlot(object = pbmc, label = TRUE) + NoLegend()
  plot_filename <- paste0("umap_", i, ".pdf")
  pdf(plot_filename)
  print(umap_plot)
  dev.off()
  cat("UMAP plot saved as:", plot_filename, "\n")
}
