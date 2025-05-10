############################
##                        ##
##  Save_atac_matrices.R  ##
##                        ##
############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--matrix',      type="character",  help='Matrix to save')
p$add_argument('--outfile',     type="character",  help='Output file')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

print(args)

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

###################
## Load metadata ##
###################

cells_metadata.dt <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

########################
## Load ArchR Project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[cells_metadata.dt$cell]

stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))

################
## PeakMatrix ##
################

if (args$matrix=="PeakMatrix") {
  
  atac.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = "PeakMatrix")
  
  # Define peak names
  row_ranges.dt <- rowRanges(atac.se) %>% as.data.table %>% 
    setnames("seqnames","chr") %>%
    .[,c("chr","start","end")] %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
  rownames(atac.se) <- row_ranges.dt$idx
  
}

#####################
## GeneScoreMatrix ##
#####################

if (grepl("GeneScoreMatrix",args$matrix)) {
  
  atac.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = args$matrix)
  
  # Define gene names
  rownames(atac.se) <- rowData(atac.se)$name
}

##########
## Save ##
##########

# Sanity checks
stopifnot(sum(duplicated(rownames(atac.se)))==0)

saveRDS(atac.se, args$outfile)