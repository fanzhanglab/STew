#' get_model_stats()
#'
#'Get Model Statistics
#'
#' This function calculates the model statistics for a given effect plot and data frame.
#'
#' @param effect_plot The effect plot containing the predicted values from different models.
#' @param dt The data frame containing the observed values.
#' @param models A character vector specifying the models to calculate statistics for.
#'
#' @return A list with two elements:
#'   - mean_var: A data frame with the mean and variance values for each gene and model.
#'   - obs_mean_var: A data frame with the observed mean and variance values for each gene.
#'
#' @export
#'
#' @import tidyverse
get_model_stats <- function(effect_plot, dt, models) {
  # Predicted
  mean_var <- effect_plot %>%
    group_by(Gene) %>%
    summarise(
      zero_mod_mean = if ("ZeroInflated" %in% models) mean(zero_mod_pred) else NA_real_,
      zero_mod_var = if ("ZeroInflated" %in% models) var(zero_mod_pred) else NA_real_,
      zero_neg_mean = if ("ZeroInflNegBinomial" %in% models) mean(zero_neg_pred) else NA_real_,
      zero_neg_var = if ("ZeroInflNegBinomial" %in% models) var(zero_neg_pred) else NA_real_,
      nb_mean = if ("NegativeBinomial" %in% models) mean(nb_pred) else NA_real_,
      nb_var = if ("NegativeBinomial" %in% models) var(nb_pred) else NA_real_,
      qua_mean = if ("Quasipoisson" %in% models) mean(qua_pred) else NA_real_,
      qua_var = if ("Quasipoisson" %in% models) var(qua_pred) else NA_real_,
      glm_mean = if ("GLM" %in% models) mean(glm_pred) else NA_real_,
      glm_var = if ("GLM" %in% models) var(glm_pred) else NA_real_,
      hurdle_mean = if ("Hurdle" %in% models) mean(hurdle_pred) else NA_real_,
      hurdle_var = if ("Hurdle" %in% models) var(hurdle_pred) else NA_real_
    )

  # Observed
  obs_mean_var =
    dt %>%
    pivot_longer(cols =  -c(nUMI
                            #,cell
    ),
    names_to = 'Gene',
    values_to = 'val') %>%
    group_by(Gene) %>%
    summarize(mean = mean(val), var = var(val))

  return(list(mean_var = mean_var, obs_mean_var = obs_mean_var))
}
