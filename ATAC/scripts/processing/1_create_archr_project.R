##############################
##                          ##
##  Create_archR_project.R  ##
##                          ##
##############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples', type="character",  nargs='+',      help='Samples')
p$add_argument('--arrow_files',     type="character",  nargs='+',      help='Arrow files')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

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