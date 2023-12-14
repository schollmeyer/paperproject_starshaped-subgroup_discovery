##################################
#  own additional functions:
##################################

gene_filter <- function(data, x_percent=6 ){
  # gene filter as described under https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/
  # Note that also genes with sdtandard deviation of epsression zero are also deleted
  i <- NULL
  for(k in seq_len(ncol(data))){
    if( mean(data[,k]>2) < x_percent/100 | mean(x[,k]>0)> 1-x_percent/100 |  sd(x[,k])==0){i <- c(i,k)}
  }
return(data[,-i])}

scaling <- function(data){
  # function that scales every column of the data matrix data by substracting the mean and dividing by the standard deviation afterwards
  for(k in seq_len(ncol(data))){
    if(sd(data[,k])==0){print("Warning: standard deviation of 0 observed")}
    data[,k] <- (data[,k]-mean(data[,k]))/sd(data[,k])
   }
return(data)}

##################################
##################################
#functions from hemberg-lab: https://github.com/hemberg-lab/scRNA.seq.datasets/blob/master/utils/create_sce.R
##################################
##################################

create_sce_from_counts <- function(counts, colData, rowData = NULL) {
    if(is.null(rowData)) {
        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                       colData = colData)
    } else {
        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                       colData = colData,
                                       rowData = rowData)
    }
    # this function writes to logcounts slot
    exprs(sceset) <- log2(calculateCPM(sceset, use.size.factors = FALSE) + 1)
    # use gene names as feature symbols
    rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if(is.null(rowData)) {
        sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
    }
    # QC
    isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    sceset <- calculateQCMetrics(sceset, feature_controls = list("ERCC" = isSpike(sceset, "ERCC")))
    return(sceset)
}

create_sce_from_normcounts <- function(normcounts, colData, rowData = NULL) {
    if(is.null(rowData)) {
        sceset <- SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)), 
                                       colData = colData)
    } else {
        sceset <- SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)), 
                                       colData = colData,
                                       rowData = rowData)
    }
    logcounts(sceset) <- log2(normcounts(sceset) + 1)
    # use gene names as feature symbols
    rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if(is.null(rowData)) {
        sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
    }
    # QC
    isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    return(sceset)
}

create_sce_from_logcounts <- function(logcounts, colData, rowData = NULL) {
    if(is.null(rowData)) {
        sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)), 
                                       colData = colData)
    } else {
        sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)), 
                                       colData = colData,
                                       rowData = rowData)
    }
    # use gene names as feature symbols
    rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if(is.null(rowData)) {
        sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
    }
    # QC
    isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    return(sceset)
}

##############################################################################
#############################################################################
