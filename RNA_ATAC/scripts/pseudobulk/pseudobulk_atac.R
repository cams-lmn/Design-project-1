#########################
##                     ##
##  Pseudobulk_atac.R  ##
##                     ##
#########################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--matrix',     type="character",         help='Matrix to pseudobulk')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--outdir',     type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE] %>% 
  dplyr::rename(cell=sample)

stopifnot(args$group_by%in%colnames(sample_metadata))

sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

table(sample_metadata[[args$group_by]])

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

print(ArchRProject)

# Sanity checks
print(args$matrix)
print(getAvailableMatrices(ArchRProject))
stopifnot(args$matrix%in%getAvailableMatrices(ArchRProject))

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(getCellColData(ArchRProject))] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
print(head(sample_metadata.to.archr))
print(str(sample_metadata.to.archr))

sample_metadata.to.archr[[args$group_by]] <- as.character(sample_metadata.to.archr[[args$group_by]])

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject))))
ArchRProject <- addCellColData(
  ArchRProject,
  data = sample_metadata.to.archr[[args$group_by]],
  name = args$group_by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

table(getCellColData(ArchRProject,args$group_by)[[1]])

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

print(sprintf("Calculating pseudobulk matrix for %s",args$matrix))

# summarise
atac.se <- getGroupSE(ArchRProject, groupBy = args$group_by, useMatrix = args$matrix, divideN = FALSE) %>% as(.,"SingleCellExperiment")

# rename features
if (grepl("peak",tolower(args$matrix),ignore.case=T)) {
  rownames(atac.se) <- rowData(atac.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
}
if (grepl("gene",tolower(args$matrix),ignore.case=T)) {
  rownames(atac.se) <- rowData(atac.se)$name
}

assayNames(atac.se) <- "counts"
head(assay(atac.se,"counts"))

assay(atac.se,"logcounts") <- log(1e6*(sweep(assay(atac.se,"counts"),2,colSums(assay(atac.se,"counts")),"/"))+1)
head(assay(atac.se,"logcounts"))

atac.se$celltype <- colnames(assay(atac.se,"counts"))
  
# save
saveRDS(atac.se, file.path(args$outdir,sprintf("%s_pseudobulk_summarized_experiment.rds",args$matrix)))