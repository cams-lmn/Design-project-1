##############################
##                          ##
##  Add_background_peaks.R  ##
##                          ##
##############################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--method',        type="character", default="chromVAR",                              help='ArchR or chromVAR')
p$add_argument('--number_background_peaks',     type="integer",    default=50,    help='Number of background peaks')
p$add_argument('--peak_metadata',     type="character",      help='Peak metadata')
p$add_argument('--genome',          type="character", default="mm10",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

print(args)

########################
## Load ArchR project ##
########################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)

#####################
## Load other data ##
#####################

args$background_peaks <- paste0(args$archr_directory,"/Background-Peaks.rds")

peak_metadata.dt <- fread(args$peak_metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

##########################
## Add background peaks ##
##########################

print("Calling background peaks!")

# This function will compute background peaks controlling for total accessibility and GC-content
# changes in the ArchR project: (1) it creates Background-Peaks.rds and (2) adds "bgdPeaks" entry to "metadata(getPeakSet(ArchRProject))"

# Background peaks are chosen by sampling peaks based on similarity in GC content and # of fragments across samples using the Mahalanobis distance. 
# - The "w" paramter controls how similar background peaks should be. 
# - The "binSize" parameter controls the precision with which the similarity is computed. Increasing "binSize" will make the function run slower.
# Returns a matrix with one row per peak and one column per iteration. values in a row represent indices of background peaks for the peak with that index

ArchRProject <- addBgdPeaks(
  ArchRProj = ArchRProject,
  nIterations = args$number_background_peaks,
  w = 0.1,
  binSize = 50,
  method = args$method,
  seed = 42,
  force = TRUE
)

#############
## Check 1 ##
#############

bgdPeaks.se <- readRDS(args$background_peaks)

tmp <- data.frame(rowRanges(bgdPeaks.se),stringsAsFactors = F)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",as.character(tmp[,1]), 
                                 as.character(tmp[,2]), as.character(tmp[,3]))

print("Check 1: before adaptation rownames. If annotation of background peaks is correct, see 'Bug Archr 1.0.2' in source code.")
print(dim(peak_metadata.dt))
print(length(rownames(bgdPeaks.se)))
print(sum(peak_metadata.dt$idx%in%rownames(bgdPeaks.se)))

#####################
## Bug Archr 1.0.2 ##
#####################

print("Adaptation backgroundpeaks for Bug Archr 1.0.2")
bgdPeaks.se <- readRDS(file.path(getOutputDirectory(ArchRProject), "Background-Peaks.rds"))
rr_bgdPeaks.se <- rowRanges(bgdPeaks.se)
df_bgdPeaks.se <- data.frame(rr_bgdPeaks.se)
print(df_bgdPeaks.se[!duplicated(df_bgdPeaks.se$seqnames),])
print(dim(df_bgdPeaks.se))
df_bgdPeaks.se$seqnames <- as.character(df_bgdPeaks.se$seqnames)
df_bgdPeaks.se$seqnames <- gsub("chr","",df_bgdPeaks.se$seqnames)
df_bgdPeaks.se_num <- df_bgdPeaks.se[!is.na(as.numeric(df_bgdPeaks.se$seqnames)),]
df_bgdPeaks.se_char <- df_bgdPeaks.se[is.na(as.numeric(df_bgdPeaks.se$seqnames)),]
df_bgdPeaks.se_num$seqnames <- as.numeric(df_bgdPeaks.se_num$seqnames)
df_bgdPeaks.se_num$seqnames <- sort(df_bgdPeaks.se_num$seqnames)
df_bgdPeaks.se <- rbind(df_bgdPeaks.se_num,df_bgdPeaks.se_char)
df_bgdPeaks.se$seqnames <- paste0(rep("chr",nrow(df_bgdPeaks.se)),df_bgdPeaks.se$seqnames)
print(df_bgdPeaks.se[!duplicated(df_bgdPeaks.se$seqnames),])
print(dim(df_bgdPeaks.se))
rr_bgdPeaks.se_new <- makeGRangesFromDataFrame(df_bgdPeaks.se,keep.extra.columns=TRUE)
rowRanges(bgdPeaks.se) <- rr_bgdPeaks.se_new
saveRDS(bgdPeaks.se,file.path(getOutputDirectory(ArchRProject), "Background-Peaks.rds"))

# or install the dev branch of ArchR don't forget to restart your R session!

#############
## Check 2 ##
#############

bgdPeaks.se <- readRDS(args$background_peaks)

tmp <- data.frame(rowRanges(bgdPeaks.se),stringsAsFactors = T)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",as.character(tmp[,1]), 
                                 as.character(tmp[,2]), as.character(tmp[,3]))

print("Check 2: after adaptation rownames. If annotation of background peaks is correct, see 'Bug Archr 1.0.2' in source code.")
print(dim(peak_metadata.dt))
print(length(rownames(bgdPeaks.se)))
print(sum(peak_metadata.dt$idx%in%rownames(bgdPeaks.se)))