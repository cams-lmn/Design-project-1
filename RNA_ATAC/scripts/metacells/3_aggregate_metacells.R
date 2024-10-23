############################
##                        ##
##  Aggregate_metacells.R ##
##                        ##
############################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata_rna',    type="character",  help='RNAseq cell metadata file')
p$add_argument('--metadata_atac',    type="character",  help='ATACseq cell metadata file')
p$add_argument('--metadata_cluster',    type="character",  help='Clustering metadata file')
p$add_argument('--archr_directory',     type="character", help='ArchR directory')
p$add_argument('--genome',          type="character", default="hg38",      help='Genome')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--sce',    type="character",  help='SingleCellExperiment')
p$add_argument('--matrices',   type="character",  nargs="+",  help='ATAC matrices')
p$add_argument('--python',   type="character",    help='Python path for reticulate')
p$add_argument('--cell2metacell',    type="character",  nargs="+", help='Metacell results')
p$add_argument('--metacell_min_reads',     type="integer",    default=1e5,     help='Minimum number of RNA reads per metacell')
p$add_argument('--metacell_min_frags',     type="integer",    default=1e5,     help='Minimum number of ATAC fragments per metacell')
p$add_argument('--incl_samples',        type="character",                               help='Suffix')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# ## START TEST ##
# args <- list()
# args$matrices <- c("PeakMatrix", "GeneScoreMatrix_TSS", "GeneScoreMatrix_distal")
# args$incl_samples <- "nodiff"
# args$metadata_rna <- "/data/homes/louisc/Project_Babraham/RNA/mapping/sample_metadata_after_doublets.txt.gz"
# args$metadata_atac <- "/data/homes/louisc/Project_Babraham/ATAC/archR/qc/sample_metadata_after_qc.txt.gz"
# args$archr_directory <- "/data/homes/louisc/Project_Babraham/ATAC/archR/"
# args$genome <- "hg38"
# args$threads <- 1
# args$sce <- "/data/homes/louisc/Project_Babraham/RNA/SingleCellExperiment.rds"
# args$matrices <- c("PeakMatrix","GeneScoreMatrix_TSS")
# args$python <- "/opt/anaconda/Anaconda3/envs/python_Babraham/bin/python" # "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
# args$metacell_min_reads <- 7.5e5
# args$metacell_min_frags <- 7.5e5
# args$filter_differentiated <- TRUE
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells"
# if (args$filter_differentiated){
#   args$cell2metacell <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/cell2metacell_assignment_nodiff.txt.gz"
# } else {
#   args$cell2metacell <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/metacells/cell2metacell_assignment.txt.gz"
# }
# ## END TEST ##

#####################
## Define settings ##
#####################

if(grepl("nodiff",args$incl_samples)){
  args$filter_differentiated <- TRUE
} else {
  args$filter_differentiated <- FALSE
}

stopifnot(file.exists(args$cell2metacell))
dir.create(args$outdir, showWarnings = F)

# Reticulate
reticulate::use_python(args$python, required = TRUE)
ad <- reticulate::import("anndata")

###########################
## Load metacell results ##
###########################

cell2metacell.dt <- args$cell2metacell %>% map(~ fread(.)) %>% rbindlist
# stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)

print(sprintf("Number of metacells: %s", length(unique(cell2metacell.dt$metacell))))

print(summary(as.numeric(table(cell2metacell.dt$metacell))))
pdf(paste0(args$outdir,"/distribution_cells_per_metacell_prefilter_",args$incl_samples,".pdf"))
hist(as.numeric(table(cell2metacell.dt$metacell)),breaks=50,main="",xlab="Amount cells per metacell")
dev.off()

###################
## Load metadata ##
###################

sample_metadata_rna <- fread(args$metadata_rna) %>%
  # .[cell%in%cell2metacell.dt$cell] %>%
  # .[,c("cell","sample","stage","genotype","celltype.mapped","closest.cell")] %>%
  # setnames("celltype.mapped","celltype") %>%
  merge(cell2metacell.dt,"cell")

sample_metadata <- data.frame(fread(args$metadata_atac)) %>%
  merge(data.frame(sample_metadata_rna))

sample_metadata <- data.table(sample_metadata)
print(head(sample_metadata))

############################################
## Filter metacells ATAC & RNAseq quality ##
############################################

nreads_per_metacell.dt <- sample_metadata[,.(nreads=sum(nCount_RNA)),by="metacell"]
nfrags_per_metacell.dt <- sample_metadata[,.(nfrags=sum(nFrags_atac)),by="metacell"]
sum(nreads_per_metacell.dt$metacell==nfrags_per_metacell.dt$metacell)
sum(nreads_per_metacell.dt$nreads>args$metacell_min_reads)
sum(nfrags_per_metacell.dt$nfrags>args$metacell_min_frags)
sum(nreads_per_metacell.dt$nreads>args$metacell_min_reads & nfrags_per_metacell.dt$nfrags>args$metacell_min_frags)

# Plot
p1 <- gghistogram(nreads_per_metacell.dt, x="nreads", y="..density..", bins=50, fill="gray70", color="gray50") +
  geom_vline(xintercept=args$metacell_min_reads, linetype="dashed") +
  labs(x="Number of reads per metacell", y="Density") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank()
  )

p2 <- gghistogram(nfrags_per_metacell.dt, x="nfrags", y="..density..", bins=50, fill="gray70", color="gray50") +
  geom_vline(xintercept=args$metacell_min_frags, linetype="dashed") +
  labs(x="Number of reads per metacell", y="Density") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank()
  )

p <- ggarrange(p1, p2)

pdf(file.path(args$outdir,paste0("nreads_metacell_threshold_",args$incl_samples,".pdf")), width=8, height=4)
print(p)
dev.off()

# Filter
print(sprintf("Removing metacells that have less than %d RNAseq reads and less than %d ATACseq fragments per metacell (%d/%d)", 
              args$metacell_min_reads, args$metacell_min_frags, 
              nrow(nreads_per_metacell.dt) - sum(nreads_per_metacell.dt$nreads>args$metacell_min_reads & nfrags_per_metacell.dt$nfrags>args$metacell_min_frags), 
              nrow(nreads_per_metacell.dt)))

metacells.to.use <- intersect(sample_metadata$metacell, 
                              nreads_per_metacell.dt$metacell[(nreads_per_metacell.dt$nreads>=args$metacell_min_reads) & 
                                                                (nfrags_per_metacell.dt$nfrags>=args$metacell_min_frags)])
sample_metadata <- sample_metadata[metacell%in%metacells.to.use]
cell2metacell.dt <- cell2metacell.dt[metacell%in%metacells.to.use]

print(summary(as.numeric(table(cell2metacell.dt$metacell))))

################################
## Filter metacells cluster 6 ##
################################

metadata_cluster <- fread(args$metadata_cluster) 

# Filter
print(sprintf("Removing metacells that are in cluster 6 (%d/%d)", 
              sum(unique(sample_metadata$metacell)%in%metadata_cluster$sample[metadata_cluster$cluster==6]),
              length(unique(sample_metadata$metacell))))

metacells.to.use <- intersect(sample_metadata$metacell, 
                              metadata_cluster$sample[metadata_cluster$cluster!=6])

sample_metadata <- sample_metadata[metacell%in%metacells.to.use]
cell2metacell.dt <- cell2metacell.dt[metacell%in%metacells.to.use]

print(table(metadata_cluster$cluster[metadata_cluster$sample%in%sample_metadata$metacell]))
print(summary(as.numeric(table(cell2metacell.dt$metacell))))

##########################################
## Aggregate RNAseq counts per metacell ##
##########################################

print(sprintf("Aggregating sequencing data for %s metacells",length(unique(sample_metadata$metacell))))

print("Calculating pseudobulk matrix for RNA")

sce <- load_SingleCellExperiment(file=args$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

sce <- sce[,colnames(sce)%in%sample_metadata$cell]

sce_pseudobulk_rna <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "metacell",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk_rna) <- "counts"

logcounts(sce_pseudobulk_rna) <- log2(1e6*(sweep(counts(sce_pseudobulk_rna),2,colSums(counts(sce_pseudobulk_rna)),"/"))+1)

##########################################
## Aggregate ATACseq frags per metacell ##
##########################################

setwd(args$archr_directory)

addArchRGenome(args$genome)
addArchRThreads(threads = args$threads) 

ArchRProject <- loadArchRProject(args$archr_directory)

# Subset metacells
ArchRProject <- ArchRProject[sample_metadata$cell]

# Sanity checks
stopifnot(args$matrices%in%getAvailableMatrices(ArchRProject))

tmp <- sample_metadata %>% 
  .[cell%in%rownames(getCellColData(ArchRProject))] %>% setkey(cell) %>% .[rownames(getCellColData(ArchRProject))] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(rownames(tmp) == rownames(getCellColData(ArchRProject))))
ArchRProject <- addCellColData(
  ArchRProject,
  data = tmp$metacell,
  name = "metacell",
  cells = rownames(tmp),
  force = TRUE
)

for (i in args$matrices) {
  print(sprintf("Calculating pseudobulk matrix for %s",i))
  
  # summarise
  atac_metacells.se <- getGroupSE(ArchRProject, groupBy = "metacell", useMatrix = i, divideN = FALSE)
  
  # rename features
  if (grepl("peak",tolower(i),ignore.case=T)) {
    rownames(atac_metacells.se) <- rowData(atac_metacells.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
  }
  if (grepl("gene",tolower(i),ignore.case=T)) {
    rownames(atac_metacells.se) <- rowData(atac_metacells.se)$name
  }
  
  # TO-DO: ADD CELL METADATA LIKE CELLTYPE
  # colData(atac_metacells.se) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame
  
  # save
  saveRDS(atac_metacells.se, file.path(args$outdir,sprintf("%s_summarized_experiment_metacells_%s.rds",i,args$incl_samples)))
}

##############################
## Define metacell metadata ##
##############################

metacell_metadata.dt <- sample_metadata %>%
  .[,c("cell","sample")] %>%
  .[cell%in%colnames(sce_pseudobulk_rna)] %>% setkey(cell) %>% .[colnames(sce_pseudobulk_rna)] %>%
  setnames("cell","metacell") %>%
  merge(cell2metacell.dt[,.(ncells=.N),by="metacell"],by="metacell")

# Add QC stats
tmp <- data.table(
  metacell = colnames(sce_pseudobulk_rna), 
  nFeature_RNA = colSums(counts(sce_pseudobulk_rna))
)
metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")

# Add QC stats
tmp <- data.table(
  metacell = colnames(atac_metacells.se), 
  nFrags_atac = colSums(assay(atac_metacells.se))
)
metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")

# Add metacell metadata as colData
stopifnot(sum(metacell_metadata.dt$metacell == colnames(sce_pseudobulk_rna))==nrow(metacell_metadata.dt))
colData(sce_pseudobulk_rna) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

saveRDS(sce_pseudobulk_rna, file.path(args$outdir,paste0("SingleCellExperiment_metacells_",args$incl_samples,".rds")))

####################################
## Plot coverage before and after ##
####################################

# to.plot <- rbind(
#   metacell_metadata.dt[,c("metacell","nFeature_RNA")] %>% .[,class:="metacell"],
#   sample_metadata[,c("metacell","nFeature_RNA")] %>% .[,class:="cell"]
# ) %>% .[,log_nFeature_RNA:=log10(nFeature_RNA)] 

# p <- gghistogram(to.plot, x="log_nFeature_RNA", fill="class", bins=100) +
#   labs(x="Number of RNA reads (log10)") +
#   theme(
#     legend.title = element_blank(),
#     axis.text = element_text(size=rel(0.75), color="black")
#   )

# pdf(file.path(args$outdir,"nreads_before_vs_after.pdf"), width=7, height=4)
# print(p)
# dev.off()

#############################################
## Convert SingleCellExperiment to anndata ##
#############################################

adata_sce <- ad$AnnData(
  X   = t(counts(sce_pseudobulk_rna)),
  obs = as.data.frame(colData(sce_pseudobulk_rna)),
  var = data.frame(gene=rownames(sce_pseudobulk_rna), row.names=rownames(sce_pseudobulk_rna))
)

# Add cell type colors
# adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
# adata_sce$uns["celltype_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(file.path(args$outdir,paste0("anndata_metacells_",args$incl_samples,".h5ad")))
fwrite(metacell_metadata.dt, file.path(args$outdir,paste0("metacells_metadata_",args$incl_samples,".txt.gz")), sep="\t", quote=F, na="NA")
