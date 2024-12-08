###########################
##                       ##
##  Create_seurat_rna.R  ##
##                       ##
###########################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--min_reads',           type="integer",    default=500,                  help='Minimum number of reads')
p$add_argument('--outdir',       type="character",                    help='Output directory')
p$add_argument('--test',          type="logical", default=FALSE,                 help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

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

count_mtx <- list()
cell.info <- list()
gene.info <- list()

for (i in args$samples) {
  print(i)
  
  
  # Load gene metadata
  gene.loc <- sprintf("%s/%s/rna_features.tsv.gz",args$inputdir,i)
  gene.info[[i]] <- fread(gene.loc, header=F)[, 1:2] %>%
    setnames(c("ens_id", "symbol")) %>%
    .[,idx:=1:.N]
  dim(gene.info[[i]])
  
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s/barcodes.tsv.gz",args$inputdir,i)
  cell.info[[i]] <- fread(barcode.loc, header=F) %>%
    setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  
  # Load matrix  
  matrix.loc <- sprintf("%s/%s/rna_matrix.mtx.gz",args$inputdir,i)
  count_mtx[[i]] <- Matrix::readMM(matrix.loc)[gene.info[[i]]$idx,]
  stopifnot(nrow(cell.info[[i]])==ncol(count_mtx[[i]]))
  rownames(count_mtx[[i]]) <- gene.info[[i]]$symbol
  colnames(count_mtx[[i]]) <- cell.info[[i]]$cell
  
}

dim_count_mtx <- lapply(count_mtx,dim)
dim_count_mtx <- do.call(cbind,dim_count_mtx)
rownames(dim_count_mtx) <- c("# features","# cells")
print(dim_count_mtx)

jpeg(paste(args$outdir,"/dim_counts_mtx.jpeg",sep=""))
plot(1:length(args$samples),dim_count_mtx[2,],type="l",xaxt="n",
     xlab="Timepoint",ylab="N cells")
axis(1,at=1:length(args$samples),labels=args$samples)
dev.off()

#######################
## Keep common genes ##
#######################

genes <- Reduce("intersect",lapply(count_mtx,rownames))
for (i in 1:length(count_mtx)) {
  count_mtx[[i]] <- count_mtx[[i]][genes,]
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
## Filter genes ##
##################

# Remove duplicated genes
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

seurat <- CreateSeuratObject(count_mtx, meta.data = cell.info.to.seurat)

head(seurat@meta.data)

# Add mitochondrial percenatge # DEBUG
mito.genes <- grep(pattern = "^mt-", x = rownames(seurat), value = TRUE, ignore.case = TRUE)
seurat[["mitochondrial_percent_RNA"]] <- PercentageFeatureSet(seurat, features = mito.genes) %>% round(2)

# Add ribosomal RNA content # DEBUG
ribo.genes <- grep(pattern = "^RP[L|S]", x = rownames(seurat), value = TRUE, ignore.case = TRUE)
seurat[["ribosomal_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA")]

##########################
## Calculate data stats ##
##########################

# foo <- apply(count_mtx[,1:5],2, function(x) table(x[x<=15])) %>% 
#   as.data.table(keep.rownames=TRUE) %>% setnames("rn","counts") %>%
#   melt(id.vars=c("counts"), variable.name="cell")

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(cell.info, paste0(args$outdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(gene.info, paste0(args$outdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))