############################
##                        ##
##  Create_arrow_files.R  ##
##                        ##
############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',           type="character",  nargs='+',      help='Samples')
p$add_argument('--fragments_files',           type="character",  nargs='+',      help='ATAC Fragments files')
p$add_argument('--genome',           type="character", default="mm10",      help='Genome')
p$add_argument('--min_fragments',     type="integer",    default=1000,   help='Minimum number of ATAC fragments')
p$add_argument('--max_fragments',     type="integer",    default=1e7,    help='Maximum number of ATAC fragments')
p$add_argument('--min_tss_score',   type="double",     default=2.5,    help='Minimum TSS score threshold')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))
print(args)

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
  logFile = paste0("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/archR/ArchRLogs/ArchR-createArrows-",args$samples,".log"),
  
  # QC metrics
  minFrags = args$min_fragments,  # The minimum number of fragments per cell
  maxFrags = args$max_fragments,  # The maximum number of fragments per cell
  minTSS = args$min_tss_score   # The minimum TSS enrichment score per cell
)