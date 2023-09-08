#' autocorrelated_genes()
#'
#' Filtering gene expression matrix by spatially auto-correlated genes and correction for multiple-testing
#'
#' @param alternative a character string specifying the alternative hypothesis, chosen from one of "greater", "two-sided" or "less". Default is "greater"
#' @param adjustPv Adjusted P value
#' @param alpha P-value threshold for LISA score to be considered significant
#' @param minPercentCells Minimum percent of cells that must be driving spatial pattern
#' @param verbose A logical parameter that controls the amount of output or information displayed during the execution of the function. When verbose = TRUE, the function will provide more detailed and informative output during its execution and when verbose = FALSE (default setting), the function will produce minimal output.
#' @param obj stCCA S3 Object
#'
#' @return Subset regularized gene expression matrix by filtered genes
#' @export
#'
#' @import MERINGUE
autocorrelated_genes <- function(obj, alternative = "greater", adjustPv = TRUE, alpha = 0.05, minPercentCells = 0.05, verbose = TRUE) {
  # Identify significantly spatial auto-correlated genes
  count <- obj$count_exp
  adjacency_matrix <- obj$adj_matrix
  I <- getSpatialPatterns(count, adjacency_matrix, alternative = alternative)

  # Filter out spatial pattern driven by small number of cells using LISA
  results.filter <- filterSpatialPatterns(mat = count,
                                          I = I,
                                          w = adjacency_matrix,
                                          adjustPv = adjustPv,
                                          alpha = alpha,
                                          minPercentCells = minPercentCells,
                                          verbose = verbose)

  # Subset regularized gene expression matrix by filtered genes
  filtered_genes <- count[results.filter,]
  obj$count_exp <- filtered_genes
  return(obj)


}
