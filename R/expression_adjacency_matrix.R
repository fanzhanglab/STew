#' To calculate cell-to-cell similarity matrix on gene (spatially auto-coorelated genes) expression based on Euclidean distance; to build up a K-nearest neighbors (KNN) graph based on the cell-to-cell similarity matrix; then to convert the KNN neighbor relation matrix into an adjacency matrix
#'
#' @param nk Number of neighbors
#' @param symm If TRUE, the connectivity matrix is made symmetrical. Default - FALSE
#' @param weight If TRUE, the weighted KNN graph is created. Default - FALSE.
#' @param mode A character string chosen from c(""directed", "undirected", "max", "min", "upper", "lower", "plus")
#' @param weighted To specify whether to create a weighted graph from an adjacency matrix. If NULL, an unweighted graph is created. If TRUE, a weighted graph is created.
#' @param obj stCCA S3 Object
#'
#' @export
#'
#'@import loe
#'Matrix
expression_adjacency_matrix <- function(obj, nk = 200, symm = FALSE, weight = TRUE, mode = "undirected", weighted = TRUE) {
  # Calculate the Euclidean distance between each pair of cells using the gene expression matrix
  mat <- t(obj$count_exp)
  distance <- dist(mat, method = "euclidean")
  distance <- as.matrix(distance)

  # Make a K-NN graph based on the Euclidean distance
  knn_matrix <- make.kNNG(distance, k = nk, symm = symm, weight = FALSE)

  # Convert K-NN neighbor relation matrix to adjacency matrix
  n <- dim(knn_matrix)[1]
  knn_adjacency <- knn_matrix

  if (!mode == "directed") {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (!weighted) {
          if ((knn_matrix[i, j] != 0 || knn_matrix[j, i] != 0)) {
            knn_adjacency[i, j] <- 1
            knn_adjacency[j, i] <- 1
          }
        } else {
          if (knn_matrix[i, j] != 0 || knn_matrix[j, i] != 0) {
            knn_adjacency[i, j] <- max(knn_matrix[i, j], knn_matrix[j, i])
            knn_adjacency[j, i] <- knn_adjacency[i, j]
          }
        }
      }
    }
  } else if (!weighted) {
    knn_adjacency[knn_adjacency != 0] <- 1
  }

  knn_adjacency <- Matrix(knn_adjacency, sparse = TRUE)
  rownames(knn_adjacency) <- rownames(mat)
  colnames(knn_adjacency) <- rownames(mat)
  obj$exp_adj_matrix <- knn_adjacency
  return(obj)
}
