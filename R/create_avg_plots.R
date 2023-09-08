#' create_avg_plots()
#'
#' Create Average Plots
#'
#' This function creates average plots based on the mean and variance data.
#'
#' @param mean_var A data frame containing the mean and variance data.
#' @param obs_mean_var A data frame containing the observed mean and variance data.
#' @param models A character vector specifying the models.
#' @param log_transform A logical value indicating whether to perform log transformation on the data.
#'
#' @return A list of average plots.
#'
#' @export
#'
#'
#' @import tidyverse
create_avg_plots <- function(mean_var, obs_mean_var, models, log_transform = FALSE) {
  if (log_transform) {
    # Log-transformed mean and variance
    plot <- obs_mean_var %>%
      ggplot(aes(x = log(mean), y = log(var))) +
      geom_point(alpha = 0.2)
  } else {
    # Mean and variance
    plot <- obs_mean_var %>%
      ggplot(aes(x = mean, y = var)) +
      geom_point(alpha = 0.2)
  }

  if ("GLM" %in% models) {
    plot <- plot +
      geom_smooth(data = mean_var %>% filter(!is.na(glm_mean)),
                  aes_string(x = if (log_transform) 'log(glm_mean)' else 'glm_mean',
                             y = if (log_transform) 'log(glm_var)' else 'glm_var',
                             color = "'GLM Poisson'"),
                  se = FALSE)
  }

  if ("LME" %in% models) {
    plot <- plot +
      geom_smooth(data = mean_var %>% filter(!is.na(lme_mean)),
                  aes_string(x = if (log_transform) 'log(lme_mean)' else 'lme_mean',
                             y = if (log_transform) 'log(lme_var)' else 'lme_var',
                             color = "'LME'"),
                  se = FALSE)
  }

  if ("ZeroInflated" %in% models) {
    plot <- plot +
      geom_smooth(data = mean_var %>% filter(!is.na(zero_mod_mean)),
                  aes_string(x = if (log_transform) 'log(zero_mod_mean)' else 'zero_mod_mean',
                             y = if (log_transform) 'log(zero_mod_var)' else 'zero_mod_var',
                             color = "'Zero Inflated'"),
                  se = FALSE)
  }

  if ("ZeroInflNegBinomial" %in% models) {
    plot <- plot +
      geom_smooth(data = mean_var %>% filter(!is.na(zero_neg_mean)),
                  aes_string(x = if (log_transform) 'log(zero_neg_mean)' else 'zero_neg_mean',
                             y = if (log_transform) 'log(zero_neg_var)' else 'zero_neg_var',
                             color = "'Zero Negative'"),
                  se = FALSE)
  }

  if ("NegativeBinomial" %in% models) {
    plot <- plot +
      geom_smooth(data = mean_var %>% filter(!is.na(nb_mean)),
                  aes_string(x = if (log_transform) 'log(nb_mean)' else 'nb_mean',
                             y = if (log_transform) 'log(nb_var)' else 'nb_var',
                             color = "'Negative Binomial'"),
                  se = FALSE)
  }

  if ("Hurdle" %in% models) {
    plot <- plot +
      geom_smooth(data = mean_var %>% filter(!is.na(hurdle_mean)),
                  aes_string(x = if (log_transform) 'log(hurdle_mean)' else 'hurdle_mean',
                             y = if (log_transform) 'log(hurdle_var)' else 'hurdle_var',
                             color = "'Hurdle'"),
                  se = FALSE)
  }

  plot <- plot +
    scale_color_discrete(name = 'Models') +
    theme_bw()

  return(plot)
}
