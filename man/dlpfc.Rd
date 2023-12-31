% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlpfc.R
\docType{data}
\name{dlpfc}
\alias{dlpfc}
\alias{DLPFC}
\title{dlpfc}
\format{
A dataframe with cell spatail coordinates and expression of highly variable genes
\describe{
\item{Spatial}{Spatial coordinates of cells}
\item{count_exp}{Gene expression of highly variable genes - a count matrix in which genes as rows with rownames and cells as columns with colnames}
}
}
\source{
\url{http://research.libd.org/spatialLIBD/}
}
\usage{
dlpfc
}
\description{
The data provided is a list of Data set containing LIBD human dorsolateral pre-frontal cortex (DLPFC)
and spatial transcriptomics data (sample id 151673) generated with the 10x Genomics Visium platform.
}
\details{
This data is provided as a sample data to understand the package and not a mandatory requirement to run the functions.
}
\keyword{datasets}
