#######################
##                   ##
##  Create_bigwig.R  ##
##                   ##
#######################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--genome',           type="character", default="mm10",      help='Genome')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--norm_method',     type="character", default="ReadsInTSS",    help='Normalisation method')
p$add_argument('--min_cells',     type="integer", default=100,    help='Minimum number of cells per celltype')
p$add_argument('--tile_size',     type="integer", default=100,    help='Tile size')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',    type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE]
colnames(sample_metadata)[1] <- "cell"

stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

# Filter celltypes by minimum number of cells
sample_metadata <- sample_metadata[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

print(table(sample_metadata[[args$group_by]]))
print(head(sample_metadata))

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

###########################
## Update ArchR metadata ##
###########################

print("addCellColData")
print(head(rownames(getCellColData(ArchRProject))))

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(getCellColData(ArchRProject))] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

print(head(sample_metadata.to.archr))

sample_metadata.to.archr$cluster <- as.character(sample_metadata.to.archr$cluster)

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

###################
## Export bigwig ##
###################

# This function will group, summarize and export a bigwig for each group in an ArchRProject.
getGroupBW(
  ArchRProj = ArchRProject,
  groupBy = args$group_by,
  normMethod = args$norm_method,
  tileSize = args$tile_size,
  maxCells = 1000, # default
  ceiling = 4
)

# Create a completion token
file.create(file.path(args$outdir,sprintf("/GroupBigWigs/%s/completed.txt",args$group_by)))
