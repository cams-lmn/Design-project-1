##################################
##                              ##
##  Dimensionality_reduction.R  ##
##                              ##
##################################

source("/kyukon/data/gent/courses/2024/design_project_C003698/groups/group01/RNA/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",  default="all",              nargs='+',     help='Samples')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed1',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--seed2',            type="integer",    default=4242,                  help='Random seed')
p$add_argument('--seed3',            type="integer",    default=424242,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output file')
p$add_argument('--vars_to_regress', type="character",   default= NULL, nargs='+',     help='Metadata columns to regress out')
p$add_argument('--batch_variable',type="character",   default="None",                            help='Metadata column to apply batch correction on')
p$add_argument("--incl_samples",  type="character",    help='Which samples should be included')
p$add_argument("--max_point_per_sample",  ttype="character",   default= NULL, help="Set a maximal value for cells per sample")
p$add_argument("--sort_samples",  type="logical",default=T, help="Should samples be sorted?")
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

dir.create(args$outdir, showWarnings = F)

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples]

# table(sample_metadata$stage)
if (args$sort_samples){
  table(sample_metadata$sample)[args$samples]
} else {
  table(sample_metadata$sample)
}
# table(sample_metadata$celltype)

#############################
## Add sequencing Run info ##
#############################

if (args$batch_variable=="SeqRun"){
  sample_metadata$SeqRun <- "Run2"
  sample_metadata$SeqRun[sample_metadata$sample=="d0"] <- "Jasmin"
  sample_metadata$SeqRun[sample_metadata$sample%in%c("d10","d14","d18","DE")] <- "Run1"
}

#####################
## Parse arguments ##
#####################

if (args$batch_variable=="None") {
  args$batch_variable <- NULL
}

#############################
## Filter naieve to primed ##
#############################

if (args$filter_differentiated){
  print("Removing differentiated cells...")
  sample_metadata <- sample_metadata[grepl("d[0-9]+",sample_metadata$sample),]
  args$samples <- args$samples[grepl("d[0-9]+",args$samples)]
  if (args$sort_samples){
    table(sample_metadata$sample)[args$samples]
  } else {
    table(sample_metadata$sample)
  }
}                                      

#################
## Filter plot ##
#################

if(!is.null(args$max_point_per_sample)){
  sample_metadata_new <- NULL
  for (j in 1:length(unique(sample_metadata$sample))){
    set.seed(args$seed3*j)
    if (sum(sample_metadata$sample==unique(sample_metadata$sample)[j])>args$max_point_per_sample){
      sample_metadata_new <- rbind(sample_metadata_new,sample_metadata[sample_metadata$sample==unique(sample_metadata$sample)[j],][sample(1:sum(sample_metadata$sample==unique(sample_metadata$sample)[j]),args$max_point_per_sample,replace=F),])
    } else {
      sample_metadata_new <- rbind(sample_metadata_new,sample_metadata[sample_metadata$sample==unique(sample_metadata$sample)[j],])
    }
  }
  sample_metadata <- sample_metadata_new
  if (args$sort_samples){
    table(sample_metadata$sample)[args$samples]
  } else {
    table(sample_metadata$sample)
  }
}

###################
## Sanity checks ##
###################

print(args$colour_by)
print(colnames(sample_metadata))
stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype) %in% names(opts$celltype.colors))

if (length(args$batch_variable)>0) {
  stopifnot(args$batch_variable%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch_variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch_variable))
    args$batch_variable <- NULL
  } else {
    library(batchelor)
  }
}

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame


#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
# if (length(args$batch_variable)>0) {
#   decomp <- modelGeneVar(sce, block=colData(sce)[[args$batch_variable]])
# } else {
#   decomp <- modelGeneVar(sce)
# }
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

if (length(args$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
  )
}

###############
## Elbowplot ##
###############

if (length(args$batch_variable)>0) {
  pca.obj <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_variable]], d = 50, 
                           preserve.single=TRUE, get.variance = TRUE)
  
  # mold pca.obj
  pca.obj.new <- pca.obj[[1]]
  colnames(pca.obj.new) <- paste(rep("PC",50),1:50,sep="")
  attr(pca.obj.new,"varExplained") <- pca.obj@metadata$var.explained
  attr(pca.obj.new,"percentVar") <- pca.obj@metadata$var.explained/pca.obj@metadata$var.total*100
  rot_obj <- pca.obj@metadata$rotation
  colnames(rot_obj) <- paste(rep("PC",50),1:50,sep="")
  attr(pca.obj.new,"rotation") <-rot_obj 
  
  reducedDim(sce_filt, "PCA") <- pca.obj.new
  jpeg(sprintf("%s/Elbowplot_batch_%s.jpg",args$outdir,args$incl_samples))
  pct_var_explained <- attr(reducedDim(sce_filt, 'PCA'), 'percentVar')
  plot(pct_var_explained)
  dev.off()
  
} else {
  sce_filt <- runPCA(sce_filt, ncomponents =50, ntop=args$features)
  jpeg(sprintf("%s/Elbowplot_%s.jpg",args$outdir,args$incl_samples))
  pct_var_explained <- attr(reducedDim(sce_filt, 'PCA'), 'percentVar')
  plot(pct_var_explained)
  dev.off()
}


############################
## PCA + Batch correction ##
############################

set.seed(args$seed1)

if (length(args$batch_variable)>0) {
  # suppressPackageStartupMessages(library(batchelor))
  # print(sprintf("Applying MNN batch correction for variable: %s", args$batch_variable))
  # outfile <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs,paste(args$batch_variable,collapse="-"))
  # pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_variable]], d = args$npcs)
  # pca.corrected <- reducedMNN(pca)$corrected
  # colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
  # reducedDim(sce_filt, "PCA") <- pca.corrected[colnames(sce),]
  
  # stopifnot(args$batch_variable=="sample")
  
  # Run PCA
  pca.obj <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_variable]], d = args$npc, preserve.single=TRUE,get.variance = TRUE)
  
  # mold pca.obj
  pca.obj.new <- pca.obj[[1]]
  colnames(pca.obj.new) <- paste(rep("PC",args$npc),1:args$npc,sep="")
  attr(pca.obj.new,"varExplained") <- pca.obj@metadata$var.explained
  attr(pca.obj.new,"percentVar") <- pca.obj@metadata$var.explained/pca.obj@metadata$var.total
  rot_obj <- pca.obj@metadata$rotation
  colnames(rot_obj) <- paste(rep("PC",args$npc),1:args$npc,sep="")
  attr(pca.obj.new,"rotation") <- rot_obj
  
  reducedDim(sce_filt, "PCA") <- pca.obj.new
  
} else {
  sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)
}

# Save PCA coordinates
pca.dt <- reducedDim(sce_filt,"PCA") %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell") # error with this 
fwrite(pca.dt, sprintf("%s/pca_features%d_pcs%d_%s.txt.gz",args$outdir, args$features, args$npcs, args$incl_samples))

##########
## UMAP ##
##########

# Run
set.seed(args$seed2)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)

# Fetch UMAP coordinates
umap.dt <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  setnames(c("UMAP1","UMAP2","cell"))

# Save UMAP coordinates
fwrite(umap.dt, sprintf("%s/umap_features%d_pcs%d_neigh%d_dist%s_%s.txt.gz",
                          args$outdir, args$features, args$npcs, args$n_neighbors, args$min_dist, args$incl_samples))

##########
## Plot ##
##########

pt.size <- ifelse(ncol(sce)>=1e4,0.8,1.2)

for (i in args$colour_by) {
  
  to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
    .[,cell:=colnames(sce_filt)] %>%
    merge(sample_metadata, by="cell")
  
  if (is.numeric(to.plot[[i]])) {
    if (max(to.plot[[i]],na.rm=T) - min(to.plot[[i]],na.rm=T) > 1000) {
      to.plot[[i]] <- log10(to.plot[[i]]+1)
      to.plot %>% setnames(i,paste0(i,"_log10")); i <- paste0(i,"_log10")
    }
  }
  
  if (args$sort_samples){
    to.plot$sample <- as.factor(to.plot$sample)
    to.plot$sample <- factor(to.plot$sample,levels=args$samples)
  }
  
  print(dim(to.plot))
  print(table(to.plot$sample))
  
  to.plot$V1 <- -to.plot$V1 
  
  p <- ggplot(to.plot, aes_string(x="V2", y="V1", fill=i)) +
    geom_point(size=pt.size, shape=21, stroke=0.05) +
    theme_classic() +
    ggplot_theme_NoAxes()
  
  # Define colormap
  if (grepl("sample",i) & !is.null(opts$color_scheme)) {
    print(cbind(args$samples,opts$color_scheme[1:length(args$samples)]))
    p <- p + scale_fill_manual(values=opts$color_scheme[1:length(args$samples)]) 
  }
  
  # Save UMAP plot
  outfile <- file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_dist%s_%s_%s.pdf",
                                           args$features, args$npcs, args$n_neighbors, args$min_dist, i, args$incl_samples))

  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
}