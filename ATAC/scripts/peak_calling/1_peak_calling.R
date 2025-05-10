#####################
##                 ##
## peak_calling.R  ##
##                 ##
#####################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--pathToMacs2',     type="character",    help='Path to MACS2 software')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--pvalue_cutoff',     type="double",   help='MACS2 p-value cutoff')
p$add_argument('--extend_summits',     type="integer",   help='Number of bp to extend peak summits')
p$add_argument('--min_cells',     type="integer",   help='Minimal amount of cells required to start peak calling')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',  type="character",                help='Output directory')
p$add_argument('--seed',     type="integer", default=42,  help='Random seed')
args <- p$parse_args(commandArgs(TRUE))

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE ]
stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

print(head(sample_metadata))

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]
ArchRProject

# # Unnecesary, but just to make sure projectMetadata is updated
tmp <- file.path(args$archr_directory,"projectMetadata.rds")
if (file.exists(tmp)) ArchRProject@projectMetadata <- readRDS(tmp)

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(getCellColData(ArchRProject))] %>% setkey(cell) %>% 
  as.data.frame() %>% tibble::column_to_rownames("cell")

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

##################
## Peak calling ##
##################

set.seed(args$seed)

ArchRProject <- addReproduciblePeakSet(
  ArchRProj = ArchRProject, 
  groupBy = args$group_by, 
  peakMethod = "Macs2",
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = args$pathToMacs2,
  cutOff = args$pvalue_cutoff,
  extendSummits = args$extend_summits,
  minCells=args$min_cells,
  plot = FALSE,
  force = TRUE
)

################
## Save peaks ##
################

print(args$outdir)

# Save PeakSet
saveRDS(ArchRProject@peakSet, file.path(args$outdir,"/PeakSet.rds"))

# fetch peaks in data.table format
dt <- getPeakSet(ArchRProject) %>% as.data.table() %>% setnames("seqnames","chr")

# Save peak metadata
fwrite(dt, file.path(args$outdir,"peak_metadata.tsv.gz"), sep="\t")

# save peaks in bed format
fwrite(dt[,c("chr","start","end")], file.path(args$outdir,"peaks_archR_macs2.bed.gz"), sep="\t", col.names = F)

#####################
## Add peak matrix ##
#####################

head(ArchRProject@peakSet)
ArchRProject@peakSet <- ArchRProject@peakSet
ArchRProject <- addPeakMatrix(ArchRProject, binarize = FALSE)
saveArchRProject(ArchRProject)