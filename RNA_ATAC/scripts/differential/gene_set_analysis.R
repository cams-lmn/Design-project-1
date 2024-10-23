###########################
##                       ##
##  Gene_set_analysis.R  ##
##                       ##
###########################

source("/data/homes/louisc/Project_Babraham/RNA_ATAC/scripts/Settings.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--DEG_overview',          type="character",    help='Differential analysis results')
p$add_argument('--outdir',          type="character",    help='Output directory')
p$add_argument('--matrix',          type="character",  nargs="+",    help='Matrix to use')
args <- p$parse_args(commandArgs(TRUE))

# args <- list()
# args$DEG_overview <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/DE_res/DEG_overview.txt"
# args$outdir <- "/data/homes/louisc/Project_Babraham/RNA_ATAC/pseudobulk/cluster/RNA/GSA"
# args$matrix <- "RNA"

####################
## Initialisation ##
####################

# Create output directory
dir.create(args$outdir)

# Set default zip commando
Sys.setenv(R_ZIPCMD="zip")

#########################################
## Load differential analysis overview ##
#########################################

DEG_results <- read.table(args$DEG_overview,sep="\t",header=T)

##############
## Universe ##
##############

# Reference genes
GenesRefFileName <- paste0(args$outdir,"/ORA_ref.txt")

write.table(DEG_results$Gene,GenesRefFileName,
            sep="\n",col.names = F,row.names = F,quote = F)

######################
## GSA marker genes ##
######################

# Genes of interest
enrichResult_all <- NULL
for (i in 1:5){
  dir_name <- paste0(args$outdir,"/Cell_Cluster",i)
  unlink(dir_name,recursive = T)
  dir.create(dir_name)
  GenesORAFileName <- paste0(dir_name,"/ORA_goi.txt")
  print(paste0("Cell Cluster ",i))
  coi <- DEG_results[,colnames(DEG_results)==paste0("clust",i,"_spec")]
  goi <- DEG_results$Gene[!is.na(coi) & coi]
  print(sprintf("%s genes of interest",length(goi)))
  if (length(goi)<50){
    print("Less than 50 differential genes so cluster is skipped for GSA")
    next()
  }
  
  write.table(goi,GenesORAFileName,
              sep="\n",col.names = F,row.names = F,quote = F)


  enrichResult_MF <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="geneontology_Molecular_Function_noRedundant",
                              interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                              referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                              dagColor = "continuous",
                              sigMethod="fdr", minNum=5, reportNum = 500,
                              projectName = "ORA_MF", outputDirectory =dir_name)

  enrichResult_BP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="geneontology_Biological_Process_noRedundant",
                              interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                              referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                              dagColor = "continuous",
                              sigMethod="fdr", minNum=5, reportNum = 500,
                              outputDirectory =dir_name, projectName = "ORA_BP")
  
  enrichResult_CC <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="geneontology_Cellular_Component_noRedundant",
                              interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                              referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                              dagColor = "continuous",
                              sigMethod="fdr", minNum=5, reportNum = 500,
                            outputDirectory =dir_name, projectName = "ORA_CC")
  
  enrichResult_KEGG <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                 enrichDatabase="pathway_KEGG",
                                 interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                                 referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                                 dagColor = "continuous",
                                 sigMethod="fdr", minNum=5, reportNum = 500,
                                 outputDirectory =dir_name, projectName = "ORA_KEGG")
  

  l_MF <- nrow(enrichResult_MF)  
  l_BP <- nrow(enrichResult_BP)  
  l_CC <- nrow(enrichResult_CC)  
  l_KEGG <- nrow(enrichResult_KEGG)  

  if (is.null(l_MF)){
    l_MF <- 0
  } 
  if (is.null(l_BP)){
    l_BP <- 0
  } 
  if (is.null(l_CC)){
    l_CC <- 0
  } 
  if (is.null(l_KEGG)){
    l_KEGG <- 0
  } 

  l_tot <- l_MF+l_BP+l_CC+l_KEGG
  print(l_MF)
  print(l_BP)
  print(l_CC)
  print(l_KEGG)
  print(l_tot)
  
  if (l_tot!=0){
    col1 <- rep(i,l_tot)
    col2 <- c(rep("MF",l_MF),rep("BP",l_BP),rep("CC",l_CC),rep("KEGG",l_KEGG))

    enrichResult_all <- rbind(enrichResult_all,
                              cbind(col1,col2,rbind(enrichResult_MF,enrichResult_BP,enrichResult_CC,enrichResult_KEGG)))
  }
}
colnames(enrichResult_all)[1:2] <- c("Cell_cluster","GOClass")
print(head(enrichResult_all))
write.table(enrichResult_all,paste0(args$outdir,"/",args$matrix,"_enrichResults_Cell_Clusters.txt"),col.names = T,row.names = F,sep="\t",quote=F)


#######################
## GSA gene_clusters ##
#######################

# Genes of interest
enrichResult_all <- NULL
gene_clusters <- unique(DEG_results$gene_cluster[!is.na(DEG_results$gene_cluster)])
for (i in gene_clusters){
  dir_name <- paste0(args$outdir,"/Gene_Cluster",i)
  unlink(dir_name,recursive = T)
  dir.create(dir_name)
  GenesORAFileName <- paste0(dir_name,"/ORA_goi.txt")
  print(paste0("Gene Cluster ",i))
  coi <- DEG_results[,colnames(DEG_results)==paste0("clust",i,"_spec")]
  goi <- DEG_results$Gene[!is.na(DEG_results$gene_cluster) & (DEG_results$gene_cluster==i)]
  print(sprintf("%s genes of interest",length(goi)))
  if (length(goi)<50){
    print("Less than 50 differential genes so cluster is skipped for GSA")
    next()
  }
  
  write.table(goi,GenesORAFileName,
              sep="\n",col.names = F,row.names = F,quote = F)


  enrichResult_MF <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="geneontology_Molecular_Function_noRedundant",
                              interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                              referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                              dagColor = "continuous",
                              sigMethod="fdr", minNum=5, reportNum = 500,
                              outputDirectory =dir_name, projectName = "ORA_MF")

  enrichResult_BP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="geneontology_Biological_Process_noRedundant",
                              interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                              referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                              dagColor = "continuous",
                              sigMethod="fdr", minNum=5, reportNum = 500,
                              outputDirectory =dir_name, projectName = "ORA_BP")
  
  enrichResult_CC <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="geneontology_Cellular_Component_noRedundant",
                              interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                              referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                              dagColor = "continuous",
                              sigMethod="fdr", minNum=5, reportNum = 500,
                            outputDirectory =dir_name, projectName = "ORA_CC")
  
  enrichResult_KEGG <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                 enrichDatabase="pathway_KEGG",
                                 interestGeneFile=GenesORAFileName, interestGeneType="genesymbol",
                                 referenceGeneFile=GenesRefFileName, referenceGeneType="genesymbol",
                                 dagColor = "continuous",
                                 sigMethod="fdr", minNum=5, reportNum = 500,
                                 outputDirectory =dir_name, projectName = "ORA_KEGG")
  
  l_MF <- nrow(enrichResult_MF)  
  l_BP <- nrow(enrichResult_BP)
  l_CC <- nrow(enrichResult_CC)
  l_KEGG <- nrow(enrichResult_KEGG)

  if (is.null(l_MF)){
    l_MF <- 0
  } 
  if (is.null(l_BP)){
    l_BP <- 0
  } 
  if (is.null(l_CC)){
    l_CC <- 0
  } 
  if (is.null(l_KEGG)){
    l_KEGG <- 0
  } 

  l_tot <- l_MF+l_BP+l_CC+l_KEGG
  print(l_MF)
  print(l_BP)
  print(l_CC)
  print(l_KEGG)
  print(l_tot)
  
  if (l_tot!=0){
    col1 <- rep(i,l_tot)
    col2 <- c(rep("MF",l_MF),rep("BP",l_BP),rep("CC",l_CC),rep("KEGG",l_KEGG))

    enrichResult_all <- rbind(enrichResult_all,
                              cbind(col1,col2,rbind(enrichResult_MF,enrichResult_BP,enrichResult_CC,enrichResult_KEGG)))
  }
}
colnames(enrichResult_all)[1:2] <- c("Gene_cluster","GOClass")
print(head(enrichResult_all))
write.table(enrichResult_all,paste0(args$outdir,"/",args$matrix,"_enrichResults_Gene_Clusters.txt"),col.names = T,row.names = F,sep="\t",quote=F)

######################
## Completion token ##
######################

# Create a completion token
file.create(paste0(args$outdir,"/completed.txt"))

