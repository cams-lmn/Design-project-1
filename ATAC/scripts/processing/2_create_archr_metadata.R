###############################
##                           ##
##  Create_archR_metadata.R  ##
##                           ##
###############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')
p$add_argument('--samples', type="character",  nargs='+',      help='Samples')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument("--sort_samples",  type="logical",default=T, help="Should samples be sorted?")
args <- p$parse_args(commandArgs(TRUE))

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)
print(dim(sample_metadata))
print(head(sample_metadata))

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

# Remove cell from metadata that were filtered out by ArchR
ArchRProject <- loadArchRProject(args$archr_directory)

print(ArchRProject)
print(rownames(ArchRProject))

# Stats unfiltered
if (args$sort_samples){
  table(ArchRProject$Sample)[args$samples]
} else {
  table(ArchRProject$Sample)
}

######################
## Load ArchR stats ##
######################

# fetch archR's metadata
archR_metadata <- getCellColData(ArchRProject) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,c("cell", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archR_metadata)%in%cols.to.rename)
colnames(archR_metadata)[idx.cols.to.rename] <- paste0(colnames(archR_metadata)[idx.cols.to.rename], "_atac")

###########
## Merge ##
###########

# print stats
print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata$cell,archR_metadata$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata$cell%in%archR_metadata$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!archR_metadata$cell%in%sample_metadata$cell)))

# merge
sample_metadata_tosave <- sample_metadata %>% 
  merge(archR_metadata,by="cell", all.y=TRUE) 

# Fill missing entries for cells that did not pass QC for the RNA
sample_metadata_tosave %>%
  .[is.na(sample),sample:=strsplit(cell,"#") %>% map_chr(1)] %>%
  .[is.na(barcode),barcode:=strsplit(cell,"#") %>% map_chr(2)]

print(head(sample_metadata_tosave))

# round
sample_metadata_tosave[,c("TSSEnrichment_atac","NucleosomeRatio_atac","PromoterRatio_atac","BlacklistRatio_atac"):=list(round(TSSEnrichment_atac,2),round(NucleosomeRatio_atac,2),round(PromoterRatio_atac,2),round(BlacklistRatio_atac,2))]
sample_metadata_tosave[,c("ribosomal_percent_RNA","mitochondrial_percent_RNA"):=list(round(ribosomal_percent_RNA,2),round(mitochondrial_percent_RNA,2))]

# sanity checks
stopifnot(all(!is.na(sample_metadata_tosave$sample)))
stopifnot(all(!is.na(sample_metadata_tosave$barcode)))

# print stats
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}

#############################
## Update ArchR's metadata ##
#############################

metadata.to.archR <- sample_metadata_tosave %>% 
  .[cell%in%rownames(getCellColData(ArchRProject))] %>% setkey(cell) %>% # .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

print(head(metadata.to.archR))

for (i in colnames(metadata.to.archR)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = metadata.to.archR[[i]], 
    name = i,
    cells = rownames(metadata.to.archR),
    force = TRUE
  )
}

print(head(getCellColData(ArchRProject)))

##########
## Save ##
##########

fwrite(sample_metadata_tosave, args$outfile, sep="\t", na="NA", quote=F)
saveArchRProject(ArchRProject, drop=TRUE)