############################
##                        ##
##  Save_atac_matrices.R  ##
##                        ##
############################

source("/data/homes/louisc/Project_Babraham/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--matrix',      type="character",  help='Matrix to save')
p$add_argument('--outfile',     type="character",  help='Output file')
p$add_argument('--genome',          type="character", default="hg38",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$archr_directory <- "/data/louisc/Project_Babraham/ATAC/archR"
# args$metadata <- "/data/louisc/Project_Babraham/ATAC/archR/qc/sample_metadata_after_qc.txt.gz"
# args$matrices <- c("PeakMatrix","GeneScoreMatrix_TSS","GeneScoreMatrix_distal")
# args$genome <- "hg38"
# args$threads <- 1
## END TEST ##

print(args)

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

###################
## Load metadata ##
###################

cells_metadata.dt <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]
# .[pass_atacQC==TRUE]

########################
## Load ArchR Project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[cells_metadata.dt$cell]

# Sanity checks
# mean(rownames(ArchRProject)%in%cells_metadata.dt$cell)
# mean(cells_metadata.dt$cell%in%rownames(ArchRProject))
# table(cells_metadata.dt[!cell%in%rownames(ArchRProject),sample])
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
  
  # Filter genes
  # atac.se <- atac.se[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm|^Mir|^Olfr",rownames(atac.se),invert=T),]
}

##########
## Save ##
##########

# Sanity checks
stopifnot(sum(duplicated(rownames(atac.se)))==0)

saveRDS(atac.se, args$outfile)
