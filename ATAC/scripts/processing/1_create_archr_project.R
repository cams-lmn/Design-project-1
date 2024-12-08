##############################
##                          ##
##  Create_archR_project.R  ##
##                          ##
##############################

source("/data/homes/louisc/Project_Babraham/ATAC/scripts/Settings.R")
suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples', type="character",  nargs='+',      help='Samples')
p$add_argument('--arrow_files',     type="character",  nargs='+',      help='Arrow files')
p$add_argument('--genome',          type="character", default="hg38",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$samples <- c("d0","d1","d3","d5","d7","d10","d14","d18","DE","NE","AME_E","AME_L")
# args$arrow_files <- paste(rep("/data/louisc/Project_Babraham/ATAC/archR/",length(args$samples)),
#                           args$samples,rep(".arrow",length(args$samples)),sep="")
# args$genome <- "hg38"
# args$threads <- 1
# args$outdir <- "/data/louisc/Project_Babraham/ATAC/archR"
## END TEST ##

#####################
## Define settings ##
#####################

# ArchR options
addArchRGenome(args$genome)

############################
## create an ArchRProject ##
############################

ArchRProject <- ArchRProject(
  ArrowFiles = args$arrow_files, 
  outputDirectory = args$outdir,
  copyArrows = FALSE
)

##########
## Save ##
##########

saveArchRProject(ArchRProject)
