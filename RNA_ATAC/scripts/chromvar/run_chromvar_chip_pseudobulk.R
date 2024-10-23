######################################
##                                  ##
##  run_chromvar_chip_pseudobulk.R  ##
##                                  ##
######################################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation')
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation')
p$add_argument('--peak_metadata',  type="character",              help='')
p$add_argument('--background_peaks',  type="character",              help='')
p$add_argument('--min_number_peaks',     type="integer",    default=30,    help='Minimum number of peaks per TF')
p$add_argument('--min_chip_score',     type="double",    default=0.15,    help='Minimum ChIP score')
p$add_argument('--ignore_negative_values',  default=TRUE,  help='Ignote negative ChIP-seq scores')
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$motif_annotation <- "CISBP"
# args$atac_peak_matrix <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/PeakMatrix_pseudobulk_summarized_experiment.rds"
# args$motifmatcher <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/virtual_chipseq/metacells/CISBP/motifmatchr_virtual_chip.rds"
# args$background_peaks <- "/data/homes/louisc/Project_Babraham/ATAC/archR/Background-Peaks.rds"
# args$peak_metadata <- "/data/homes/louisc/Project_Babraham/ATAC/archR/peak_metadata.tsv.gz"
# args$min_number_peaks <- 30
# args$min_chip_score <- 0.30
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/chromvar_chip/pseudobulk"
# args$ignore_negative_values <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

dir.create(args$outdir, showWarnings=F, recursive=T)

if (grepl("pseudobulk",args$outdir)){
  args$type <- "pseudobulk"
} else {
  args$type <- "cells"
}

##########################
## Load ATAC PeakMatrix ##
##########################

atac_peakMatrix.se <- readRDS(args$atac_peak_matrix)#[,opts$celltypes]

# Load peak metadata
peak_metadata.dt <- fread(args$peak_metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("idx","score","GC")] %>% 
  setkey(idx) %>% .[rownames(atac_peakMatrix.se)]
stopifnot(peak_metadata.dt$idx==rownames(atac_peakMatrix.se))

print(atac_peakMatrix.se)
print(head(peak_metadata.dt))

# temporary
if (any(!c("start","strand")%in%colnames(rowData(atac_peakMatrix.se)))) {
  tmp <- rownames(atac_peakMatrix.se) %>% strsplit(":") %>% map_chr(2)
  rowData(atac_peakMatrix.se)$start <- tmp %>% strsplit("-") %>% map_chr(1)
  rowData(atac_peakMatrix.se)$end <- tmp %>% strsplit("-") %>% map_chr(2)
}

###############################
## Load motifmatcher results ##
###############################

print("Loading motifmatcher results...")

motifmatcher.se <- readRDS(args$motifmatcher)

stopifnot(c("motifMatches","VirtualChipScores")%in%names(assays(motifmatcher.se)))

###################################################################
## Update motifmatcher results using the virtual ChIP-seq scores ##
###################################################################

if (args$ignore_negative_values){
  print(sprintf("Number of matches before filtering negative TF binding values: %d",sum(assay(motifmatcher.se,"motifMatches"))))
  # assays(motifmatcher.se) <- assays(motifmatcher.se)[""]
  assay(motifmatcher.se,"motifMatches")[assay(motifmatcher.se,"VirtualChipScores")<0] <- F
  print(sprintf("Number of matches after filtering negative TF binding values: %d",sum(assay(motifmatcher.se,"motifMatches"))))
}

print(sprintf("Number of matches before filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher.se,"motifMatches"))))
# assays(motifmatcher.se) <- assays(motifmatcher.se)[""]
assay(motifmatcher.se,"motifMatches")[abs(assay(motifmatcher.se,"VirtualChipScores"))<=args$min_chip_score] <- F
print(sprintf("Number of matches after filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher.se,"motifMatches"))))

assays(motifmatcher.se) <- assays(motifmatcher.se)["motifMatches"]

################
## Filter TFs ##
################

# Filter TFs with too few peaks
TFs <- which(colSums(assay(motifmatcher.se,"motifMatches"))>=args$min_number_peaks) %>% names
TFs.removed <- which(colSums(assay(motifmatcher.se,"motifMatches"))<args$min_number_peaks) %>% names

cat(sprintf("%s TFs removed because they don't have enough binding sites: %s", length(TFs.removed), paste(TFs.removed, collapse=" ")))

motifmatcher.se <- motifmatcher.se[,TFs]
motifmatcher.se
head(assay(motifmatcher.se,"motifMatches")[,1:6])


###########################
## Load background peaks ##
###########################

print("Loading background peaks...")

bgdPeaks.se <- readRDS(args$background_peaks)
tmp <- rowRanges(bgdPeaks.se)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
stopifnot(sort(rownames(bgdPeaks.se))==sort(rownames(atac_peakMatrix.se)))
bgdPeaks.se <- bgdPeaks.se[rownames(atac_peakMatrix.se),]
bgdPeaks.se
head(assay(bgdPeaks.se,"bgdPeaks")[,1:6])


###################
## Sanity checks ##
###################

nrow(motifmatcher.se)
nrow(atac_peakMatrix.se)
sum(rownames(atac_peakMatrix.se)==rownames(motifmatcher.se))
sum(rownames(atac_peakMatrix.se)%in%rownames(motifmatcher.se))
stopifnot(sort(rownames(atac_peakMatrix.se))==sort(rownames(motifmatcher.se)))

#########################################
## Use chromVAR default implementation ##
#########################################

print("Running chromVAR's default implementation...")

# prepare data for chromvar
assayNames(atac_peakMatrix.se) <- "counts"
summary(assay(atac_peakMatrix.se))

# Compute deviations
chromvar_deviations_chromvar.se <- computeDeviations(
  object = atac_peakMatrix.se,
  annotations = motifmatcher.se,
  background_peaks = assay(bgdPeaks.se)
)
chromvar_deviations_chromvar.se 
summary(assay(chromvar_deviations_chromvar.se,"z"))

# Save
saveRDS(chromvar_deviations_chromvar.se, file.path(args$outdir,sprintf("chromVAR_deviations_chip_%s_%s_chromvar.rds",
                                                                       args$motif_annotation,args$type)))


chromvar_variability_chromvar.se <- computeVariability(
  object = chromvar_deviations_chromvar.se,
)

pdf(file.path(args$outdir,sprintf("chromVAR_variability_chip_%s_%s_chromvar.pdf",args$motif_annotation,args$type)))
plotVariability(chromvar_variability_chromvar.se, use_plotly = FALSE) 
dev.off()

#########################################
## Use ArchR's chromVAR implementation ##
#########################################

print("Running ArchR's chromVAR implementation...")

# source(here::here("atac/archR/chromvar/utils.R"))

featureDF <- data.frame(
  rowSums = rowSums(assay(atac_peakMatrix.se)),
  start = rowData(atac_peakMatrix.se)$start,
  end = rowData(atac_peakMatrix.se)$end,
  GC = peak_metadata.dt$GC
)

chromvar_deviations_archr.se <- .customDeviations(
  countsMatrix = assay(atac_peakMatrix.se),
  annotationsMatrix = as(assay(motifmatcher.se),"dgCMatrix"),
  backgroudPeaks = assay(bgdPeaks.se),
  expectation = featureDF$rowSums/sum(featureDF$rowSums),
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
)
chromvar_deviations_archr.se 
summary(assay(chromvar_deviations_archr.se,"z"))

# Save
saveRDS(chromvar_deviations_archr.se, file.path(args$outdir,sprintf("chromVAR_deviations_chip_%s_%s_archr.rds",
                                                                    args$motif_annotation,args$type)))
