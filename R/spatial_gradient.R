#' spatial_gradient()
#'
#' Spatial-informed cell gradient plots
#'
#' @param t A character string to define plot title. Null in default.
#' @param obj stCCA S3 Object
#'
#' @export
#'
#' @import ggplot2
#' ggpubr
spatial_gradient <- function(obj, t = NULL) {

  coordis <- obj$spatial
  joint_embedding <- obj$joint_emb

  if(all(rownames(coordis) != rownames(joint_embedding))) stop("Spatial coordinates should have the same cell barcode as canonical variate")
  coordis2 <- cbind(joint_embedding[,1:20],coordis)

  coordi_plots <- list()
  for (i in colnames(coordis2)[1:20]) {

    coordi_plots[[i]] <- ggplot(coordis2, aes(x= coordis[, 2],y= coordis[, 1])) +
      geom_point(aes(color= .data[[i]])) + #+geom_point(shape = 1,colour = "black")
      scale_color_gradient2(low = "#4dac26", mid = "#f0f0f0", high = "#d01c8b", midpoint=0) +
      theme_bw() +
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
      labs(title = t) +
      theme(plot.title = element_text(size=14,face="bold"))  + scale_y_reverse()


  }

  return(coordi_plots)
}
