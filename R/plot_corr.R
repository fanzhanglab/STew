#' Barplot of canonical correlation for each pair of canonical variate
#'
#' @param cc Output of CCA() function. Result of Sparse CCA.
#' @param t Title
#' @param tsub Subtitle
#'
#' @export
#'
#' @import ggplot2
#' ggpubr
plot_corr <- function(cc, t, tsub=NULL) {

  cca_cor <- data.frame(cc$cors)
  cca_cor$CC <- paste0(1:20, sep="")
  colnames(cca_cor) <- c("cor", "CC")

  bar <- ggbarplot(cca_cor, x = "CC", y = "cor", fill = "#beaed4", color = "#beaed4", ylim = c(0,1), xlab = "Dimension", ylab = "Canonical correlation coefficients", cex.names=.5, main = t) +
    labs(subtitle = tsub)
  return(bar)

}
