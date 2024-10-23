##########################
##                      ##
##  Settings.R - Joint  ##
##                      ##
##########################


#####################
# Loading libraries #
#####################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(gganimate))
suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(xfun))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(Rclusterpp))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(motifmatchr))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(chromVAR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(JASPAR2020))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(WebGestaltR))

#############
# Functions #
#############

# load_SingleCellExperiment
############

load_SingleCellExperiment <- function(file, normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE) {
  library(SingleCellExperiment); library(scran); library(scater);
  sce <- readRDS(file)
  if (!is.null(cells)) sce <- sce[,cells]
  if (!is.null(features)) sce <- sce[features,]
  if (remove_non_expressed_genes) sce <- sce[which(Matrix::rowSums(counts(sce))>15),]
  if (normalise) sce <- logNormCounts(sce)
  
  return(sce)
}

# ggplot_theme_NoAxes
############

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

# tfidf
############

tfidf <- function(mtx, method = 1, scale.factor = 1e4) {
  npeaks <- colSums(mtx)
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }
  
  tf <- Matrix::tcrossprod(mtx, y = Matrix::Diagonal(x=1/npeaks))
  
  rsums <- rowSums(mtx)
  if (any(rsums == 0)) {
    warning("Some features contain 0 total counts")
  }
  idf <- ncol(mtx) / rsums
  
  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    tf <- log1p(tf * scale.factor)
    idf <- log(1 + idf)
  }
  mtx.tfidf <- Matrix::Diagonal(n = length(idf), x = idf) %*% tf
  
  if (method == 1) {
    mtx.tfidf <- log1p(mtx.tfidf * scale.factor)
  }
  colnames(mtx.tfidf) <- colnames(mtx)
  rownames(mtx.tfidf) <- rownames(mtx)
  
  # set NA values to 0
  mtx.tfidf[is.na(mtx.tfidf)] <- 0
  
  return(mtx.tfidf)
}


matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

.augment_matrix <-function(mtx, samples) {
  samples <- unique(samples)
  mtx <- t(mtx)
  aug_mtx<-matrix(NA, ncol=ncol(mtx), nrow=length(samples))
  aug_mtx<-mtx[match(samples,rownames(mtx)),,drop=FALSE]
  rownames(aug_mtx)<-samples
  colnames(aug_mtx)<-colnames(mtx)
  return(t(aug_mtx))
}

minmax.normalisation <- function(x)
{
  min_x <- min(x,na.rm=T)
  max_x <- max(x,na.rm=T)
  res <- (x-min_x)/(max_x-min_x)
  print("Minmax report")
  print(sprintf("All values lie between %s and %s.",min(x,na.rm=T),max(x,na.rm=T)))
  print(sprintf("Selected boundaries for minmax normalisation are %s and %s.",min_x,max_x))
  return(res)
}

minmax.10pct.normalisation <- function(x)
{
  min_x <- quantile(x,0.05,na.rm=T)
  max_x <- quantile(x,0.95,na.rm=T)

  res <- (x-min_x)/(max_x-min_x)
  res[res<0] <- 0
  res[res>1] <- 1
  print("Minmax report")
  print(sprintf("All values lie between %s and %s.",min(x,na.rm=T),max(x,na.rm=T)))
  print(sprintf("Selected boundaries for minmax normalisation are %s and %s.",min_x,max_x))
  return(res)
}

minmax.fb.normalisation <- function(x,lb,ub)
{
  res <- (x-lb)/(ub-lb)
  res[res<0] <- 0
  res[res>1] <- 1
  return(res)
}

getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
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

# z_score
############

z_score <- function(x){
  (x - mean(x)) / sd(x)
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

# Utils chromVAR
############

.customDeviations <- function(
  countsMatrix = NULL,
  annotationsMatrix = NULL,
  backgroudPeaks = NULL,
  expectation = NULL,
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
  ){

  # sanity checks
  stopifnot(nrow(countsMatrix) == nrow(backgroudPeaks))
  stopifnot(length(expectation) == nrow(countsMatrix))

  colData <- DataFrame(seq_len(ncol(countsMatrix)), row.names = colnames(countsMatrix))[,FALSE]
  norm_expectation <- expectation / sum(expectation) #Double check this sums to 1!
  countsPerSample <- Matrix::colSums(countsMatrix)

  d <- max(floor(ncol(annotationsMatrix)/20), 1)
  m <- 0
  results <- .safelapply(seq_len(ncol(annotationsMatrix)), function(x){
    if(x %% d == 0){
      m <- 1 #Print to console
    }
    if(x %% max(floor(d/5), 2) == 0){
      if(m != 1){
      }else{
        m <- 0 #Reset
      }
    }
    if(x %% max(c(d, 10)) == 0){
      gc()
    }
    .customDeviationsSingle(
      annotationsVector = annotationsMatrix[, x, drop=FALSE],
      countsMatrix = countsMatrix,
      backgroudPeaks = backgroudPeaks,
      countsPerSample = countsPerSample,
      expectation = norm_expectation,
      out = out,
      prefix = prefix
    )
  }, threads = threads)
  cn <- colnames(countsMatrix)
  rm(countsMatrix)
  gc()

  # parse output
  if("z" %in% tolower(out)){
    z <- t(vapply(results, function(x) x[["z"]], rep(0, length(cn))))
  }else{
    z <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  if("deviations" %in% tolower(out)){
    dev <- t(vapply(results, function(x) x[["dev"]], rep(0, length(cn))))
  }else{
    dev <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  colnames(z) <- cn
  colnames(dev) <- cn

  #Check First
  nullOverlap <- is.null(results[[1]]$overlap)
  rowData <- lapply(seq_along(results), function(x){
      resx <- results[[x]]
      if(nullOverlap){
        data.frame(fractionMatches = resx$matches)
      }else{
        data.frame(fractionMatches = resx$matches, fractionBackgroundOverlap = resx$overlap)
      }
    }) %>% Reduce("rbind",.)
  rownames(rowData) <- colnames(annotationsMatrix)

  # craete output summarized experiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = dev, z = z),
    colData = colData,
    rowData = rowData
  )
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)[tolower(out)]
  return(se)
}

.customDeviationsSingle <- function(
  annotationsVector = NULL,
  countsMatrix = NULL,
  countsPerSample = NULL,
  backgroudPeaks = NULL,
  out = c("deviations", "z"),
  expectation = NULL,
  threshold = 1,
  prefix = ""
  ){
  .binarizeMat <- function(mat = NULL){
    mat@x[mat@x > 0] <- 1
    mat
  }

  if (length(annotationsVector@x) == 0) {
      out <- list(
        z = rep(NA, ncol(countsMatrix)),
        dev = rep(NA, ncol(countsMatrix)),
        expFG = NA,
        expBG = NA,
        matches = 0,
        overlap = NA
        )
    return(out)
  }

  outList <- tryCatch({

    ################################
    # Fore Ground Deviations
    ################################
    observed <- as.vector(Matrix::t(annotationsVector) %*% countsMatrix)
    expected <- as.vector(Matrix::t(annotationsVector) %*% expectation %*% countsPerSample)
    observed_deviation <- (observed - expected)/expected

    #Filter those with no matches at all
    fail_filter <- which(expected == 0)

    ################################
    # Back Ground Deviations
    ################################
    if("z" %in% tolower(out)){

      #Compute Background Null Per Iteration
      niterations <- ncol(backgroudPeaks)
      sampleMat <- Matrix::sparseMatrix(
          j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
          i = rep(seq_len(niterations), each = length(annotationsVector@x)),
          x = rep(annotationsVector@x, niterations),
          dims = c(niterations, nrow(countsMatrix))
      )
      sampled <- as.matrix(sampleMat %*% countsMatrix)
      sampledExpected <- sampleMat %*% expectation %*% countsPerSample
      sampledDeviation <- (sampled - sampledExpected)/sampledExpected
      bgOverlap <- Matrix::mean(.binarizeMat(sampleMat) %*% .binarizeMat(annotationsVector)) / length(annotationsVector@x)

      #Summary
      meanSampledDeviation <- Matrix::colMeans(sampledDeviation)
      sdSampledDeviation <- apply(as.matrix(sampledDeviation), 2, sd)

      #Norm Deviation
      normdev <- (observed_deviation - meanSampledDeviation)
      z <- normdev/sdSampledDeviation
      if (length(fail_filter) > 0) {
        z[fail_filter] <- NA
        normdev[fail_filter] <- NA
      }

    }else{

      #Compute Background Null Per Iteration
      niterations <- ncol(backgroudPeaks)
      sampleMat2 <- Matrix::sparseMatrix(
          j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
          i = rep(1, niterations * length(annotationsVector@x)),
          x = rep(annotationsVector@x, niterations),
          dims = c(1, nrow(countsMatrix))
      )
      sampled2 <- (sampleMat2 %*% countsMatrix)[1,]
      sampledExpected2 <- (sampleMat2 %*% expectation %*% countsPerSample)[1,]
      ######################
      # Equivalent to above
      # colMeans(sampled) - colMeans(sampledExpected))/colMeans(sampledExpected)
      ######################
      sampledDeviation2 <- (sampled2 - sampledExpected2)/sampledExpected2
      bgOverlap <- NA

      #Norm Deviation
      normdev <- (observed_deviation - sampledDeviation2)
      z <- NULL
      if (length(fail_filter) > 0) {
        normdev[fail_filter] <- NA
      }

    }

    outList <- list(
      z = z,
      dev = normdev,
      matches = length(annotationsVector@x) / nrow(countsMatrix),
      overlap = bgOverlap
    )
    outList
  }, error = function(e){
    errorList <- list(
      annotationsVector = annotationsVector,
      observed = if(exists("observed", inherits = FALSE)) observed else "observed",
      expected = if(exists("expected", inherits = FALSE)) expected else "expected",
      sampleMat = if(exists("sampleMat", inherits = FALSE)) sampleMat else "sampleMat",
      sampleMat2 = if(exists("sampleMat", inherits = FALSE)) sampleMat2 else "sampleMat2",
      sampledDeviation = if(exists("sampledDeviation", inherits = FALSE)) sampledDeviation else "sampledDeviation",
      sampledDeviation2 = if(exists("sampledDeviation2", inherits = FALSE)) sampledDeviation2 else "sampledDeviation2",
      normdev = if(exists("normdev", inherits = FALSE)) normdev else "normdev",
      z = if(exists("z", inherits = FALSE)) z else "z"
    )
  })
  return(outList)
}


.safelapply <- function(..., threads = 1, preschedule = FALSE){
  if(threads > 1){
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    errorMsg <- list()

    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }

    if(length(errorMsg) != 0){
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
    }
  } else{
    o <- lapply(...)
  }
  o
}

# give.n
############

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

# theme_graph
############

# theme_graph <- function() {
#   theme_classic() +
#
#     theme(
#       base_family = "Arial Narrow" ,
#       axis.line = element_blank(),
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       axis.title = element_blank(),
#       legend.position = "none"
#     )
#
# }


# output_plot
############

output_plot <- function(p,fn,width,height,UMAP=FALSE){
  axes_line_size <- 1.5

  pdf(paste0(fn,".pdf"), width=width, height=height)
  print(p)
  dev.off()

  png(paste0(fn,".png"), width=width, height=height, units="in", res=300)
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

  png(paste0(fn,"_nolabs.png"), width=width, height=height, units="in", res=300)
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
# opts$celltype.colors <- c("#DAF7A6","#FFC300","#DDAA05","#C70039","#900C3F","#EBECF0")

opts$naieve_pluri <- c("KLF4","KLF5","TFCP2L1","DNMT3L","FGF4","KLF17")
opts$general_pluri <- c("NANOG","POU5F1","SALL4","SOX2","TDGF1")
opts$common_post_impl <- c("ETV4","ETV5","MYC","SOX11","FZD7","CDH2","SALL2","SFRP2","ZIC2")
opts$housekeeping <- c("GADPH","ACTB","TBP","B2M","HPRT1","HMBS","PGK1","PPIA","GUSB","TFRC","YWHAZ","SDHA","LDHA","UBC")


