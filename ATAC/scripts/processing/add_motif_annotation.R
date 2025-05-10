##############################
##                          ##
##  Add_motif_annotation.R  ##
##                          ##
##############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',     type="character", help='ArchR directory')
p$add_argument('--peak_calls',     type="character", help='Peak calls file (.rds)')
p$add_argument('--cutoff',     type="double",    default=1e-4,    help='Motif match p-value cutoff')
p$add_argument('--width',     type="integer",    default=7,    help='Minimum motif length')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--motif_annotation',     type="character", default="JASPAR", help='Motif annotation (CISBP or JASPAR)')
p$add_argument('--folder_manual_motifs',     type="character", default=NULL, help='Folder with manually defined motifs encoded in a .jaspar format')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################
  
# I/O
io <- list()
io$positions_outfile <- file.path(args$archr_directory, sprintf("Annotations/%s-Positions.rds",args$motif_annotation))   
io$scores_outfile <- file.path(args$archr_directory, sprintf("Annotations/%s-Scores.rds",args$motif_annotation))
io$motif2gene_outfile <- file.path(args$archr_directory, sprintf("Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
io$archR.projectMetadata <- file.path(args$archr_directory,"projectMetadata.rds")
io$peakAnnotations.output <- file.path(args$archr_directory, "Annotations/peakAnnotation.rds")

########################
## Load ArchR Project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads) 

ArchRProject <- loadArchRProject(args$archr_directory)

# Load ArchR projectMetadata
ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)

# Load peaks
ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(args$peak_calls), force = TRUE)
names(ArchRProject@peakSet) <- sprintf("%s:%s-%s",seqnames(ArchRProject@peakSet), start(ArchRProject@peakSet), end(ArchRProject@peakSet))

# Sanity checks
ArchRProject

#########################
## Load motif database ##
#########################

print("Loading motif annotations...")

if (grepl("JASPAR",toupper(args$motif_annotation))) {
  library(JASPAR2020)
  
  # this returns a PFMatrixList
  # Mus musculus: 107 TF motifs
  motifs <- getMatrixSet(JASPAR2020, list(species = "Mus musculus", collection = "CORE"))
  
  # this returns a list with two entries: (1) motifs (PFMatrixList), (2) motifSummary (a data.frame). 
  # Note TFs are renamed to capitalised gene symbol (FOXF2_1, FOXD1_2, IRF2_3)
  obj <- .summarizeJASPARMotifs(motifs)
  motifs <- obj$motifs
  motifSummary <- obj$motifSummary
  motifSummary$name <- NULL
  
  # Filter out TF motifs
  tmp <- rownames(obj$motifSummary)
  tmp <- tmp[grep("\\.\\.",tmp,invert=TRUE)] # fusion proteins
  motifs <- motifs[tmp]
  motifSummary <- motifSummary[tmp,]
  
  # Rename TF motifs
  motifSummary$symbol <- strsplit(rownames(motifSummary),"\\_") %>% map_chr(1)
  rename <- c("(NKX[0-9])\\.([0-9])" = "\\1-\\2"); motifSummary$symbol <- stringr::str_replace_all(motifSummary$symbol,rename)
  motifSummary$symbol <- strsplit(motifSummary$symbol,"\\.") %>% map_chr(1)
  
  stopifnot(sum(duplicated(motifSummary$name))==0)
  
} else if (grepl("CISBP",toupper(args$motif_annotation))) {
  library(chromVARmotifs)
  
  data("mouse_pwms_v2")
  motifs <- mouse_pwms_v2 
  obj <- .summarizeChromVARMotifs(motifs)
  motifs <- obj$motifs
  motifSummary <- obj$motifSummary
  
  # Filter out TF motifs
  
  # Rename TF motifs
  motifSummary$symbol <- motifSummary$name; motifSummary$name <- NULL
} else {
  stop("Motif annotation not recognised")
}

stopifnot(names(motifs)==rownames(motifSummary))

########################
## Add manual matches ##
########################

# Save motif2gene
motif2gene.dt <- motifSummary %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>% .[,strand:=NULL] %>% setnames("symbol","gene")
fwrite(motif2gene.dt, io$motif2gene_outfile, sep="\t", quote=F)

####################################
## Run motifmatchr: get positions ##
####################################

print("Running motifmatchr to get positions...")

motifmatcher_positions.out <- matchMotifs(
  pwms = motifs,
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  out = "positions",
  p.cutoff = args$cutoff,
  w = args$width
)

# Save
saveRDS(motifmatcher_positions.out, io$positions_outfile)

#################################
## Run motifmatchr: get scores ##
#################################

print("Running motifmatchr to get scores...")

motifmatcher_scores.out <- matchMotifs(
  pwms = motifs,
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  out = "scores",
  p.cutoff = args$cutoff,
  w = args$width
)

# print
print("MotifScores:"); assay(motifmatcher_scores.out,"motifScores")[1:5,1:5]
print("MotifMatches:"); assay(motifmatcher_scores.out,"motifMatches")[1:5,1:5]
print("MotifCounts:"); assay(motifmatcher_scores.out,"motifCounts")[1:5,1:5]
print(sprintf("Matches with cutoff=%s and width=%s: %.2f%%", args$cutoff, args$width, 100*mean(assay(motifmatcher_scores.out,"motifMatches")>0)))

# Save
saveRDS(motifmatcher_scores.out, io$scores_outfile)

###########################
## Save motif annotation ##
###########################

print("Saving peakAnnotations file...")

if (file.exists(io$peakAnnotations.output)) {
  peakAnnotation <- readRDS(io$peakAnnotations.output)
} else {
  # peakAnnotation <- list()
  peakAnnotation <- S4Vectors::SimpleList()
  
}

if (args$motif_annotation%in%names(peakAnnotation)) {
  print(sprintf("%s already exists in peakAnnotation.rds, replacing...",args$motif_annotation))
}

peakAnnotation[[args$motif_annotation]]$Name <- args$motif_annotation
peakAnnotation[[args$motif_annotation]]$motifs <- motifs
peakAnnotation[[args$motif_annotation]]$motifSummary <- motifSummary
peakAnnotation[[args$motif_annotation]]$Positions <- io$positions_outfile
peakAnnotation[[args$motif_annotation]]$Matches <- io$scores_outfile

saveRDS(peakAnnotation, io$peakAnnotations.output)
