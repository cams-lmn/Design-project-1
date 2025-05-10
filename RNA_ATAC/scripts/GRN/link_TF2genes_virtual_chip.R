####################################
##                                ##
##  Link_TF2genes_virtual_chip.R  ##
##                                ##
####################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation')
p$add_argument('--virtual_chip_matrix',       type="character",                help='TF2peak virtual chip matrix')
p$add_argument('--peak2gene',       type="character",                help='Links between peaks and genes based on genomic distance')
p$add_argument('--max_distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimum in silico ChIP-seq score')
p$add_argument('--outfile',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###########################
## Load virtual ChIP-seq ##
###########################

virtual_chip.mtx <- readRDS(args$virtual_chip_matrix)
print(head(virtual_chip.mtx))

#########################################################
## Load peak2gene linkages using only genomic distance ##
#########################################################

peak2gene.dt <- fread(args$peak2gene) %>%
  .[dist<=args$max_distance] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

# Sanity checks
stopifnot(length(intersect(rownames(virtual_chip.mtx),unique(peak2gene.dt$peak)))>1e4)

#########################################################
## Link TFs to target genes using the virtual ChIP-seq ##
#########################################################

tf2gene_chip.dt <- colnames(virtual_chip.mtx) %>% map(function(i) {
  print(i)
  # Select target peaks (note that we only take positive correlations into account)
  target_peaks_i <- names(which(assay(virtual_chip.mtx,"VirtualChipScores")[,i]>=args$min_chip_score))
  
  if (length(target_peaks_i)>=1) {
    
    tmp <- data.table(
      tf = i,
      peak = target_peaks_i,
      chip_score =  assay(virtual_chip.mtx,"VirtualChipScores")[target_peaks_i,i]
    ) %>% merge(peak2gene.dt[peak%in%target_peaks_i,c("peak","gene","dist")], by="peak")
    
    return(tmp)
  }
}) %>% rbindlist

# Save
fwrite(tf2gene_chip.dt, args$outfile, sep="\t")
