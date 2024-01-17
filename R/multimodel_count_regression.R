#' multimodel_count_regression()
#'
#' Multimodel Count Regression
#'
#' This function performs multimodel count regression analysis.
#'
#' @param outcomes_interest A character vector specifying the outcomes of interest.
#' @param models A character vector specifying the models to be fitted.
#' @param plot_density A logical value indicating whether to plot density plots of the fitted models.
#' @param plot_qq A logical value indicating whether to plot Q-Q plots of the fitted models.
#' @param estimate A logical value indicating whether to estimate model parameters and standard errors.
#' @param avg_effect A logical value indicating whether to compute average marginal effects.
#' @param rootogram A logical value indicating whether to generate rootograms for the fitted models.
#' @param full_model The formula for the full model.
#' @param null_model The formula for the null model.
#' @param data The data frame containing the variables of interest.
#' @param workers The number of parallel workers to use for parallel computing.
#' @param path The path to save the plots.
#'
#' @return A list containing the results of the multimodel count regression analysis.
#'
#' @export
#' @import future
#' future.apply
#' princurve
#' tidyverse
#' data.table
#' parallel
#' MASS
#' lme4
#' lmtest
#' future
#' pscl
#' vcd
#' ComplexHeatmap
#' ggpubr
#' gridExtra
#' ggrepel
#' viridis
#' parameters
multimodel_count_regression <- function(outcomes_interest,  # a vector of genes of interest
                                        models,          # a vector of models to use (e.g., "GLM", "ZeroInflated", etc.)
                                        plot_density,    # a logical value, if TRUE, plots density
                                        plot_qq,         # a logical value, if TRUE, plots Q-Q plot
                                        estimate,        # a logical value, if TRUE, calculates estimates
                                        avg_effect,        # a logical value, if TRUE, plots the effect of the expression over the genes.
                                        rootogram,       # a logical value, if TRUE, creates rootogram
                                        full_model,      # a formula describing the full model to be tested
                                        null_model,      # a formula describing the null model to be compared with the full model
                                        data,      # a submatrix of the data set to use for the analysis, MUST BE COUNT DATA NOT LOG NORMALIZED
                                        workers,
                                        path) {       # number of workers for parallel computation

  # Set up the parallel environment
  plan(multisession, workers = workers)

  # Apply the main function to each element of genes_interest in parallel
  results <- future_lapply(outcomes_interest,
                           FUN = function(x) count_model_main(x,
                                                              models,
                                                              plot_density,
                                                              plot_qq,
                                                              estimate,
                                                              log_plot = avg_effect,
                                                              rootogram,
                                                              full_model,
                                                              null_model,
                                                              sub_matrix = data,
                                                              path))
  return(results)
}
