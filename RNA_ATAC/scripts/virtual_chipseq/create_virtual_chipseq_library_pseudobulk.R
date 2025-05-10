###################################################
##                                               ##
##  Create_virtual_chipseq_library_pseudobulk.R  ##
##                                               ##
###################################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--tf2peak_cor',  type="character",              help='Correlations between TF RNA expression and peak accessibility')
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)')
p$add_argument('--peak_metadata',             type="character",        help='Peak metadata file')
p$add_argument('--motif2gene',  type="character",              help='Motif annotation')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation')
p$add_argument('--min_number_peaks',     type="integer",    default=50,    help='Minimum number of peaks per TF')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

# Parallel processing
# plan(sequential)
plan(multisession, workers=args$threads)

###############################
## Load pseudobulk ATAC data ##
###############################

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)
atac_peakMatrix_pseudobulk.se 
rowData(atac_peakMatrix_pseudobulk.se) <- NULL
atac_peakMatrix_pseudobulk.se 

# Load peak metadata
peak_metadata.dt <- fread(args$peak_metadata) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

######################
## Load motifmatchR ##
######################

motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_peakMatrix_pseudobulk.se)))
motifmatcher.se <- motifmatcher.se[rownames(atac_peakMatrix_pseudobulk.se),]

###############################
## Load TF2peak correlations ##
###############################

tf2peak_cor.se <- readRDS(args$tf2peak_cor)

stopifnot(rownames(tf2peak_cor.se)%in%rownames(motifmatcher.se))
stopifnot(rownames(tf2peak_cor.se)%in%rownames(atac_peakMatrix_pseudobulk.se))

##################
## Subset peaks ##
##################

peaks <- Reduce(intersect,list(rownames(motifmatcher.se),rownames(tf2peak_cor.se),peak_metadata.dt$idx))

print(sprintf("Number of peaks: %s",length(peaks)))

tf2peak_cor.se <- tf2peak_cor.se[peaks,]
atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[peaks,]
peak_metadata.dt <- peak_metadata.dt[idx%in%peaks]

################################
## Load motif2gene annotation ##
################################

motif2gene.dt <- fread(args$motif2gene)

################
## Rename TFs ##
################

print("Rename TFs")

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
TFs <- intersect(colnames(tf2peak_cor.se),motif2gene.dt$gene)
length(motifs)
length(TFs)

motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%TFs]
motifs <- motif2gene_filt.dt$motif
TFs <- motif2gene_filt.dt$gene
length(motifs)
length(TFs)

tmp <- TFs
names(tmp) <- motifs

stopifnot(motif2gene_filt.dt$motif%in%colnames(motifmatcher.se))
stopifnot(motif2gene_filt.dt$gene%in%colnames(tf2peak_cor.se))

motifmatcher.se <- motifmatcher.se[,motifs]
colnames(motifmatcher.se) <- tmp[colnames(motifmatcher.se)]
tf2peak_cor.se <- tf2peak_cor.se[,TFs]
stopifnot(colnames(motifmatcher.se)==colnames(tf2peak_cor.se))

##################
## Prepare data ##
##################

tf2peak_cor.mtx <- assay(tf2peak_cor.se,"cor")
motifmatcher.mtx <- assay(motifmatcher.se,"motifScores")
atac.mtx <- assay(atac_peakMatrix_pseudobulk.se[peaks,],"logcounts") %>% round(3)

######################################
## Create virtual chip-seq library ##
######################################

print("Predicting TF binding sites...")

opts$TFs <- TFs

stopifnot(!duplicated(opts$TFs))
print(sprintf("Number of TFs: %s",length(opts$TFs)))

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]

virtual_chip.dt <- opts$TFs %>% map(function(i) {
  print(sprintf("%s (%s/%s)",i,match(i,opts$TFs),length(opts$TFs)))
  
  peaks <- names(which(abs(tf2peak_cor.mtx[,i])>0)) 
  
  if (length(peaks)>=args$min_number_peaks) {
    
    # calculate accessibility score
    max_accessibility_score <- apply(atac.mtx[peaks,],1,max) %>% round(2)
    
    # calculate correlation score
    correlation_score <- tf2peak_cor.mtx[peaks,i] %>% round(2)
    correlation_score[correlation_score==0] <- NA
    
    # calculate motif score
    motif_score <- motifmatcher.mtx[peaks,i]
    motif_score <- round(motif_score/max(motif_score),2)
    
    predicted_score <- correlation_score * minmax.normalisation(max_accessibility_score*motif_score)
    
    tmp <- data.table(
      peak = peaks, 
      correlation_score = correlation_score,
      max_accessibility_score = max_accessibility_score,
      motif_score = motif_score,
      score = round(predicted_score,2)
    ) %>% sort.abs("score") %>%
      .[,c("peak","score","correlation_score","max_accessibility_score","motif_score")]
    
    bed.dt <- tmp %>% copy %>%
      .[,peak:=str_replace(peak,":","-")] %>%
      .[,chr:=strsplit(peak,"-") %>% map_chr(1)] %>%
      .[,start:=strsplit(peak,"-") %>% map_chr(2)] %>%
      .[,end:=strsplit(peak,"-") %>% map_chr(3)] %>%
      .[,c("chr","start","end","score")]
    
    # Save 
    data_dir <- file.path(args$outdir, "data")
    if (!dir.exists(data_dir)) {
      dir.create(data_dir, recursive = TRUE)
    }
    fwrite(tmp, sprintf("%s/data/%s.txt.gz",args$outdir,i), sep="\t", quote=F, col.names = T)
    fwrite(bed.dt, sprintf("%s/data/%s.bed.gz",args$outdir,i), sep="\t", quote=F, col.names = F)
    
    to_return.dt <- tmp[!is.na(score),c("peak","score")] %>% .[,tf:=i]
    return(to_return.dt)
  }
}) %>% rbindlist

########################
## Save sparse matrix ##
########################

virtual_chip.mtx <- virtual_chip.dt %>% 
  .[,peak:=factor(peak,levels=rownames(motifmatcher.se))] %>%
  data.table::dcast(peak~tf, value.var="score", fill=0, drop=F) %>%
  matrix.please %>% Matrix::Matrix(.)

saveRDS(virtual_chip.mtx, file.path(args$outdir,"virtual_chip_matrix.rds"))

###################################################################
## Update motifmatchr results using the virtual ChIP-seq library ##
###################################################################

print("Updating motifmatchr results using the virtual ChIP-seq library...")

motifmatcher_chip.se <- motifmatcher.se[,colnames(virtual_chip.mtx)]

# reset motif matches
stopifnot(rownames(virtual_chip.mtx)==rownames(motifmatcher_chip.se))
stopifnot(colnames(virtual_chip.mtx)==colnames(motifmatcher_chip.se))
assay(motifmatcher_chip.se,"VirtualChipScores") <- virtual_chip.mtx
saveRDS(motifmatcher_chip.se, file.path(args$outdir,"motifmatchr_virtual_chip.rds"))