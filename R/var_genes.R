#' Calculating row mean for each gene using normalized gene expression matrix
#'
#' @param obj stCCA object
#'
#' @export
#'
#'@import tibble
#'Matrix
#'matrixStats
var_genes <- function(obj){

  environment(FindVariableGenesSeurat) <- asNamespace("Seurat", base.OK = TRUE)

  exp <- Matrix(obj$count_exp, sparse = TRUE)
  vargenes <- FindVariableGenesSeurat(exp)
  vargenes_means_sds <- tibble(symbol = rownames(exp), gene_info = vargenes, mean = Matrix::rowMeans(exp))
  vargenes_means_sds$stddev <- matrixStats::rowSds(as.matrix(exp))
  obj$hvg <- vargenes_means_sds

  return(obj)

}

