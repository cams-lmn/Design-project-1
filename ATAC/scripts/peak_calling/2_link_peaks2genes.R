##########################
##                      ##
##  link_peaks2genes.R  ##
##                      ##
##########################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/ATAC/scripts/Settings.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--gene_metadata',  type="character",                help='Gene metadata')
p$add_argument('--peak_metadata',  type="character",                help='Peak metadata')
p$add_argument('--gene_window',  type="integer",  default=1e5,               help='Genomic window size')
p$add_argument('--outdir',  type="character",                help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F)

###############
## Load data ##
###############

# Load gene metadata
gene_metadata <- fread(args$gene_metadata) %>%
    setnames(old = c("Gene stable ID", "Gene name", "Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)", "Strand"), 
             new = c("ens_id", "gene", "chr", "start", "end", "strand")) %>%
    .[, chr := as.factor(sub("chr", "", chr))] %>%
    .[, strand := ifelse(as.numeric(strand) == 1, "+", "-")] %>%
    .[, c("chr", "start", "end", "gene", "ens_id", "strand")]

print(head(gene_metadata))

# Load peak metadata
peakSet.dt <- fread(args$peak_metadata) %>%
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,c("chr","start","end")] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,start,end)] %>%
  setkey(chr,start,end)

print(head(peakSet.dt))

#############
## Overlap ##
#############

gene_metadata.ov <- copy(gene_metadata) %>%
  .[strand=="+",c("gene.start","gene.end"):=list(start,end)] %>%
  .[strand=="-",c("gene.start","gene.end"):=list(end,start)] %>%
  .[strand=="+",c("start","end"):=list (gene.start-args$gene_window, gene.end+args$gene_window)] %>%
  .[strand=="-",c("end","start"):=list (gene.start+args$gene_window, gene.end-args$gene_window)] %>% 
  setkey(chr,start,end)

stopifnot((gene_metadata.ov$end-gene_metadata.ov$start)>0)

ov <- foverlaps(
  peakSet.dt,
  gene_metadata.ov,
  nomatch = NA
) %>%  .[,c("start","end"):=NULL] %>%
  setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
  .[,peak.mean:=(peak.start+peak.end)/2] %>%
  # calculate distance from the peak to the genebody
  .[,dist:=min(abs(gene.end-peak.mean), abs(gene.start-peak.mean)), by=c("gene","ens_id","peak","strand")] %>%
  .[strand=="+" & peak.mean>gene.start & peak.mean<gene.end,dist:=0] %>%
  .[strand=="-" & peak.mean<gene.start & peak.mean>gene.end,dist:=0]

# Select nearest gene
ov_nearest <- ov %>%
  .[.[,.I[dist==min(dist)], by=c("peak")]$V1] %>%
  .[complete.cases(.)] %>%
  .[!duplicated(peak)]

##########
## Save ##
##########

fwrite(ov, file.path(args$outdir,"peaks2genes_all.txt.gz"), sep="\t", na="NA")
fwrite(ov_nearest, file.path(args$outdir,"peaks2genes_nearest.txt.gz"), sep="\t", na="NA")