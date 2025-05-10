##############################
##                          ##
##  add_GeneScore_matrix.R  ##
##                          ##
##############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)

#####################
## Define gene set ##
#####################

genes.gr <- getGenes(ArchRProject)

#########################################
## Add Gene Scores using default model ##
#########################################

print("Adding GeneScores distal")

# Note that this will add the matrices to the arrowFiles
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  geneModel = "exp(-abs(x)/5000) + exp(-1)", # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_distal",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),  # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 5000,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 0,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 500,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 5,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM"),
  force = TRUE
)

#############################################
## Add Gene Scores ignoring distal regions ##
#############################################

print("Adding GeneScores TSS")

# TSS
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  useTSS = TRUE,
  extendTSS = TRUE,
  geneModel = "1",                    # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_TSS",
  extendUpstream = c(0,0),
  extendDownstream = c(0,0),          # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 500,                 # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 100,               # Number of bp downstream the gene to extend the gene body.
  tileSize = 100,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 1,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM"),
  force = TRUE
)

# Create a completion token
file.create("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/archR/addGeneScoreMatrix_completed.txt")