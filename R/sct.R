#' sct()
#'
#' Apply variance stabilizing transformation to UMI count data using a regularized Negative Binomial regression model. This will remove unwanted effects from UMI data
#'
#' @param n_genes Number of genes
#' @param obj stCCA S3 Object
#'
#' @export
#'
#' @import sctransform
#'
sct <- function(obj, n_genes = 2000) {
  umi <- obj$count_exp
  original_colnames <- colnames(umi)
  original_rownames <- rownames(umi)
  umi = umi %>% as.matrix()
  colnames(umi) = original_colnames
  rownames(umi) = original_rownames

  sct_count = sctransform::vst(umi, n_genes = n_genes, verbosity = 1)$y
  sct_count <- as.matrix(sct_count)
  obj$count_exp <- sct_count
  return(obj)

}
