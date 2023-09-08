#' plot_cluster()
#'
#' Visualizing clusters in human brain image
#'
#' @param coordis Spatial location of cells with cell barcode as rows and coordinates as columns.
#' @param label Cluster label as a column of data frame.
#' @param colors Color palette for clusters.
#' @param t Plot Title. Null in default.
#' @param subt Plot subtitle. Null in default.
#' @export
#'
#' @import ggplot2
#' ggpubr
plot_cluster <- function(coordis, label, colors, t = NULL, subt = NULL) {

  coordis$label <- label


  p <- coordis %>%
    ggplot() +
    geom_point(aes(x = coordis[, 2], y = coordis[, 1], color = factor(label))) +
    scale_color_manual(values = colors, name = '') +
    theme(axis.line = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank())

}

