############################
##                        ##
##  Settings.R - ATACseq  ##
##                        ##
############################


####################
## Load libraries ##
####################

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(gganimate))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Rclusterpp))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(chromVAR))
suppressPackageStartupMessages(library(motifmatchr))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))


##############
# Functions #
#############

# ggplot_theme_NoAxes
#############

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}


# summarize motif functions
############

.summarizeJASPARMotifs <- function(motifs = NULL){
  
  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)
  
  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
      family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
      alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  
  names(motifs) <- motifNames
  
  out <- list(motifs = motifs, motifSummary = motifDF)
  
  return(out)
  
}

.summarizeChromVARMotifs <- function(motifs = NULL){
  
  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)
  
  motifNames2 <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex
  }) %>% unlist(.)
  
  motifDF <- lapply(seq_along(motifs), function(x){
    df <- data.frame(
      row.names = motifNames[x],
      name = motifNames2[[x]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  
  names(motifs) <- motifNames
  
  out <- list(motifs = motifs, motifSummary = motifDF)
  
  return(out)
  
}

# dropNA
############

dropNA <- function(x) {
  if(!is(x, "matrix")) stop("x needs to be a matrix!")
  
  zeros <- which(x==0, arr.ind=TRUE)
  ## keep zeros
  x[is.na(x)] <- 0
  x[zeros] <- NA
  x <- Matrix::drop0(x)
  x[zeros] <- 0
  x
}


# Pseudobulk
############

pseudobulk_sce_fn <- function(x, assay = NULL, by, fun = c("sum", "mean", "median"), scale = FALSE) {

  # check validity of input arguments
  fun <- match.arg(fun)
  if (is.null(assay))  assay <- assayNames(x)[1]

  # store aggregation parameters &
  # nb. of cells that went into aggregation
  md <- metadata(x)
  md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)

  # get aggregation function
  # fun <- switch(fun, sum = "rowSums", mean = "rowMeans", median = "rowMedians")

  # drop missing factor levels & tabulate number of cells
  cd <- dplyr::mutate_if(as.data.frame(colData(x)), is.factor, droplevels)
  colData(x) <- DataFrame(cd, row.names = colnames(x),check.names = FALSE)
  md$n_cells <- table(as.data.frame(colData(x)[, by]))

  # assure 'by' colData columns are factors so that missing combinations aren't dropped
  for (i in by)
    if (!is.factor(x[[i]]))
      x[[i]] <- factor(x[[i]])

  # split cells & compute pseudo-bulks
  cs <- .split_cells(x, by)
  # pb <- .pb(x, cs, assay, fun)
  pb <- .pb(x=x, by=by, fun=fun)
  if (scale & length(by) == 2) {
    ls <- lapply(.pb(x, cs, "counts", "rowSums"), colSums)
    pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
    names(pb) <- names(ls)
  }

  # construct SCE
  pb <- SingleCellExperiment(pb, metadata = md)

  # propagate 'colData' columns that are unique across 2nd 'by'
  if (length(by) == 2) {
    cd <- colData(x)
    ids <- colnames(pb)
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(x))))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  return(pb)
}


.split_cells <- function(x, by) {
  if (is(x, "SingleCellExperiment")) x <- colData(x)
  cd <- data.frame(x[by], check.names = FALSE)
  cd <- data.table(cd, cell = rownames(x)) %>% split(by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(cd, length(by), "cell")
}

.pb <- function(x, by, fun) {

  # compute pseudobulks
  # y <- scuttle::summarizeAssayByGroup(x, assay.type = assay, ids = (ids <- colData(x)[by]), statistics = fun, BPPARAM = BiocParallel::SerialParam())
  y <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[by], statistics = fun)
  colnames(y) <- y[[by[length(by)]]]

  if (length(by) == 1)  return(assay(y))

  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
    ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- map(is, ~assay(y)[, .])

  # fill in missing combinations
  for (i in seq_along(ys)) {
    fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
    if (length(fill != 0)) {
      foo <- matrix(0, nrow(x), length(fill))
      colnames(foo) <- fill
      foo <- cbind(ys[[i]], foo)
      o <- paste(sort(unique(y[[by[2]]])))
      ys[[i]] <- foo[, o]
    }
  }
  return(ys)
}


# output_plot
############

output_plot <- function(p,fn,width,height,UMAP=FALSE){
  axes_line_size <- 1.5

  pdf(paste0(fn,".pdf"), width=width, height=height)
  print(p)
  dev.off()
  
  if (UMAP){
    p_nolabs <- p +
      theme(
        legend.position = "none",
        plot.title=element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x =  element_blank(),
        axis.line=element_blank()
      ) 
      width <- 7
      height <- 7
  } else {
    p_nolabs <- p +
      theme(
        legend.position = "none",
        plot.title=element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x =  element_blank(),
        axis.line=element_line(size=axes_line_size),
        axis.ticks=element_line(size=axes_line_size,colour="black"),
        axis.ticks.length=unit(.25, "cm"),
      )  
  }


  pdf(paste0(fn,"_nolabs.pdf"), width=width, height=height)
  print(p_nolabs)
  dev.off()
}

###############
# Set options #
###############

opts <- list()
opts$color_scheme <- c("#A06932","#DBB216","#EFE32A","#D7DB54","#A7B019","#9AA126",
                       "#7E7721","#5D6821","#353D8A","#CA3639","#DCADCD","#7A327E")

opts$celltype.colors <- c("#DAF7A6","#FFC300","#BCA043","#C70039","#900C3F","#EBECF0","#FFFFFF")