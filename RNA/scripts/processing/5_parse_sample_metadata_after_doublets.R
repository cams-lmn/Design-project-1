##############################################
##                                          ##
##  parse_sample_metadata_after_doublets.R  ##
##                                          ##
##############################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',         type="character",   help='Metadata file to use as input')
p$add_argument('--doublet_files',    type="character", nargs="+",  help='Results of the doublet score detection algorithm')
p$add_argument('--outfile',          type="character",   help='Output file')
args <- p$parse_args(commandArgs(TRUE))

##########################
## Load mapping results ##
##########################

doublet.dt <- args$doublet_files %>% map(~ fread(.)) %>% rbindlist

# stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)
# stopifnot(mapping_seurat.dt$cell%in%sample_metadata$cell)

####################
## Merge and save ##
####################

to.save <- fread(args$metadata) %>% 
  merge(doublet.dt[,c("cell","scDblFinder.score","doublet_call")] %>% setnames("scDblFinder.score","doublet_score"), by="cell", all.x=TRUE)
fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)

# to.save[pass_rnaQC==T & is.na(celltype.mapped_mnn)]
# stopifnot(to.save[pass_rnaQC==T & is.na(celltype.mapped_mnn),.N]==0)
# stopifnot(to.save[pass_rnaQC==T & is.na(celltype.mapped_seurat),.N]==0)