############################
##                        ##
##  Create_arrow_files.R  ##
##                        ##
############################

source("/data/homes/louisc/Project_Babraham/ATAC/scripts/Settings.R")
suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',           type="character",  nargs='+',      help='Samples')
p$add_argument('--fragments_files',           type="character",  nargs='+',      help='ATAC Fragments files')
p$add_argument('--genome',           type="character", default="hg38",      help='Genome')
p$add_argument('--min_fragments',     type="integer",    default=1000,   help='Minimum number of ATAC fragments')
p$add_argument('--max_fragments',     type="integer",    default=1e7,    help='Maximum number of ATAC fragments')
p$add_argument('--min_tss_score',   type="double",     default=2.5,    help='Minimum TSS score threshold')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))
print(args)

## START TEST ##
# args <- list()
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$fragments_files <- paste(rep("/data/louisc/Project_Babraham/results_sierra/",length(args$samples)),
#                               args$samples,rep("/atac_fragments.tsv.gz",length(args$samples)),sep="")
# args$genome <- "hg38"
# args$min_fragments <- 1000
# args$max_fragments <- 1e7
# args$min_tss_score <- 2.5
# args$threads <- 1
# args$outdir <- "/data/louisc/Project_Babraham/ATAC/archR"
## END TEST ##

#####################
## Define settings ##
#####################

setwd(args$outdir)

# ArchR options
addArchRThreads(threads=args$threads) 
addArchRGenome(args$genome)

rhdf5::h5disableFileLocking()

########################
## create Arrow Files ##
########################

# did not work for all files together?
# ArrowFiles <- rep("",length(args$samples))
# for (i in 1:length(args$samples)){
#   ArrowFiles[i] <- createArrowFiles(
#     inputFiles = args$fragments_files[i],
#     sampleNames = args$samples[i],
#     outputNames = args$samples[i],
#     addTileMat = FALSE,
#     addGeneScoreMat = FALSE,
#     excludeChr = c("chrM", "chrY"),
#     
#     subThreading = FALSE, # parallel processing doesn't work well (https://github.com/GreenleafLab/ArchR/issues/248)
#     force = TRUE,
#     
#     # QC metrics
#     minFrags = args$min_fragments,  # The minimum number of fragments per cell
#     maxFrags = args$max_fragments,  # The maximum number of fragments per cell
#     minTSS = args$min_tss_score   # The minimum TSS enrichment score per cell
#   )
# }

ArrowFiles <- createArrowFiles(
  inputFiles = args$fragments_files,
  sampleNames = args$samples,
  outputNames = args$samples,
  addTileMat = FALSE,
  addGeneScoreMat = FALSE,
  excludeChr = c("chrM", "chrY"),
  
  subThreading = FALSE, # parallel processing doesn't work well (https://github.com/GreenleafLab/ArchR/issues/248)
  force = TRUE,
  cleanTmp = FALSE,
  logFile = paste0("/data/homes/louisc/Project_Babraham/ATAC/archR/ArchRLogs/ArchR-createArrows-",args$samples,".log"),
  
  # QC metrics
  minFrags = args$min_fragments,  # The minimum number of fragments per cell
  maxFrags = args$max_fragments,  # The maximum number of fragments per cell
  minTSS = args$min_tss_score   # The minimum TSS enrichment score per cell
)

warnings()