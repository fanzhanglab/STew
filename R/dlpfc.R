#' dlpfc
#'
#'
#'  The data provided is a list of Data set containing LIBD human dorsolateral pre-frontal cortex (DLPFC)
#'  and spatial transcriptomics data (sample id 151673) generated with the 10x Genomics Visium platform.
#'
#' This data is provided as a sample data to understand the package and not a mandatory requirement to run the functions.
#'
#'@docType data
#' @name dlpfc
#' @aliases DLPFC
#'
#' @source \url{http://research.libd.org/spatialLIBD/}
#' @format A dataframe with cell spatail coordinates and expression of highly variable genes
#' \describe{
#'   \item{Spatial}{Spatial coordinates of cells}
#'   \item{count_exp}{Gene expression of highly variable genes - a count matrix in which genes as rows with rownames and cells as columns with colnames}
#' }

"dlpfc"
