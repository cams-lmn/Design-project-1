################################
##                            ##
## archR_add_GroupCoverage.R  ##
##                            ##
################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--min_cells',     type="integer",    default=50,   help='Minimum number of cells')
p$add_argument('--max_cells',     type="integer",    default=1000,   help='Maximum number of cells')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]
stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

# Filter celltypes by minimum number of cells
sample_metadata <- sample_metadata[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

print(head(sample_metadata))

########################
## Load ArchR project ##
########################

print(args$genome)

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

print(ArchRProject)
print(head(rownames(getCellColData(ArchRProject))))

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(getCellColData(ArchRProject))] %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject))))
ArchRProject <- addCellColData(
  ArchRProject,
  data = sample_metadata.to.archr[[args$group_by]],
  name = args$group_by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

# print cell numbers
table(getCellColData(ArchRProject,args$group_by)[[1]])

#########################
## Add Group Coverages ##
#########################

# This function will merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
# Output: creates files in archR/GroupCoverages/celltype: [X]._.Rep[Y].insertions.coverage.h5
ArchRProject <- addGroupCoverages(ArchRProject, 
                                  groupBy = args$group_by,
                                  useLabels = FALSE,  # do not use sample information
                                  minCells = args$min_cells,
                                  maxCells = args$max_cells,
                                  force = TRUE
)

##########
## Save ##
##########

saveRDS(ArchRProject@projectMetadata, paste0(args$archr_directory,"/projectMetadata.rds"))