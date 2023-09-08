library(future)
library(future.apply)
library(princurve)
library(tidyverse)
library(Seurat)
library(data.table)
library(parallel)
library(MASS)
library(lme4)
library(lmtest)
library(future)
library(pscl)
library(vcd)
library(ComplexHeatmap)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(viridis)
library(parameters)
count_model_main <- function(x,
                            models,
                            plot_density = FALSE,
                            plot_qq = FALSE,
                            estimate = FALSE,
                            log_plot = FALSE,
                            rootogram = FALSE,
                            full_model,
                            null_model,
                            sub_matrix,
                            path = getwd()){

  meta_colors <- list(
    "model" = c(
      "Bayes Negative Binomial" = '#A6CEE3',
      "Bayes Poisson" = '#1F78B4',
      "Negative Binomial" = '#B2DF8A',
      "Poisson Intercept Mixed Effect" = '#33A02C',
      "Zero Inflated Negative Binomial" = '#FB9A99',
      "Zero Inflated Poisson" = '#E31A1C',
      'Hurdle' = 'blue'
    )
  )

  ## Model run

  results_list <- list()

  if ("GLM" %in% models){
    tryCatch({
      mod = glm(formula(paste0(x, full_model)), data = sub_matrix, family = poisson)
      mod1 = glm(formula(paste0(x, null_model)), data = sub_matrix, family = poisson)
    }, error = function(e) {
      NULL
    })
  }

  if ("LME" %in% models){
    tryCatch({
      lme_mod = lme4::glmer(formula(paste0(x, full_model)), data = sub_matrix , family = 'poisson')
      lme_mod1 = lme4::glmer(formula(paste0(x, null_model)), data = sub_matrix , family = 'poisson')
    }, error = function(e) {
      NULL
    })
  }

  if ("ZeroInflated" %in% models){
    tryCatch({
      zero_mod = zeroinfl(formula(paste0(x, full_model)), data = sub_matrix, dist = 'poisson')
      zero_mod1 = zeroinfl(formula(paste0(x, null_model)), data = sub_matrix , dist = 'poisson')
    }, error = function(e) {
      NULL
    })
  }

  if ("ZeroInflNegBinomial" %in% models){
    tryCatch({
      zero_neg = zeroinfl(formula(paste0(x, full_model)), data = sub_matrix , dist = 'negbin')
      zero_neg1 = zeroinfl(formula(paste0(x, null_model)), data = sub_matrix , dist = 'negbin')
    }, error = function(e) {
      NULL
    })
  }

  if ("NegativeBinomial" %in% models){
    tryCatch({
      nb = glm.nb(formula(paste0(x, full_model)), data = sub_matrix)
      nb1 = glm.nb(formula(paste0(x, null_model)), data = sub_matrix)
    }, error = function(e) {
      NULL
    })
  }

  if ("Hurdle" %in% models){
    tryCatch({
      hurdle_mod =  hurdle(formula(paste0(x, full_model)), data = sub_matrix, dist = "poisson", zero.dist = "binomial")
      hurdle_mod1 =  hurdle(formula(paste0(x, null_model)), data = sub_matrix, dist = "poisson", zero.dist = "binomial")
    }, error = function(e) {
      NULL
    })
  }

  if ("Quasipoisson" %in% models){
    tryCatch({
      qua = glm(formula(paste0(x, full_model)), data = sub_matrix, family = 'quasipoisson')
      qua1 = glm(formula(paste0(x, null_model)), data = sub_matrix, family = 'quasipoisson')
    }, error = function(e) {
      NULL
    })
  }

  ## plot density

  if (plot_density == TRUE) {
    gg <- ggplot()

    if ("GLM" %in% models) {
      gg <- gg +
        geom_density(aes(x = predict(mod, type = 'response'), color = 'Poisson Intercept Mixed Effect'), alpha = 0.2)
    }

    if ("LME" %in% models) {
      gg <- gg +
        geom_density(aes(x = predict(lme_mod, type = 'response'), color = 'Poisson Intercept Mixed Effect'), alpha = 0.2)
    }

    if ("ZeroInflated" %in% models) {
      gg <- gg +
        geom_density(aes(x = predict(zero_mod, type = 'response'), color = 'Zero Inflated Poisson'), alpha = 0.2)
    }

    if ("ZeroInflNegBinomial" %in% models) {
      gg <- gg +
        geom_density(aes(x = predict(zero_neg, type = 'response'), color = 'Zero Inflated Negative Binomial'), alpha = 0.2)
    }

    if ("NegativeBinomial" %in% models) {
      gg <- gg +
        geom_density(aes(x = predict(nb, type = 'response'), color = 'Negative Binomial'), alpha = 0.2)
    }

    if ("Hurdle" %in% models) {
      gg <- gg +
        geom_density(aes(x = predict(hurdle_mod, type = 'response'), color = 'Hurdle'), alpha = 0.2)
    }

    gg <- gg +
      geom_histogram(aes(x = unlist(sub_matrix[x]), y = ..density.., fill = 'Observed'), alpha = 0.2) +
      scale_color_discrete(name = 'Models') +
      xlab('log(mRNA count)') +
      scale_fill_grey(name = '') +
      ggtitle(x) +
      theme_bw() +
      ylab('Density') +
      scale_x_continuous(trans = 'log10') +
      xlab('log(CPM + 1)') +
      scale_color_manual(values = meta_colors$model)


    gg
    setwd(path)
    ggsave(paste0(x, '_model_comparison', Sys.Date(), '.png'))
  }


  ## qq plot
  if (plot_qq == TRUE) {
    df <- data.frame()
    tp <- list()
    sub <- character()

    models_names <- c("GLM", "LME", "ZeroInflated", "ZeroInflNegBinomial", "NegativeBinomial", "Hurdle")
    models_labels <- c("Poisson Intercept Mixed Effect", "Zero Inflated Poisson", "Zero Inflated Negative Binomial", "Negative Binomial", "Hurdle")

    if ("GLM" %in% models) {
      tryCatch({
        tp$GLM <- residuals(mod)
      }, error = function(e) {
        tp$GLM <- NULL
      })
    }

    if ("LME" %in% models) {
      tryCatch({
        tp$LME <- residuals(lme_mod)
      }, error = function(e) {
        tp$LME <- NULL
      })
    }

    if ("ZeroInflated" %in% models) {
      tryCatch({
        tp$ZeroInflated <- residuals(zero_mod)
      }, error = function(e) {
        tp$ZeroInflated <- NULL
      })
    }

    if ("ZeroInflNegBinomial" %in% models) {
      tryCatch({
        tp$ZeroInflNegBinomial <- residuals(zero_neg)
      }, error = function(e) {
        tp$ZeroInflNegBinomial <- NULL
      })
    }

    if ("NegativeBinomial" %in% models) {
      tryCatch({
        tp$NegativeBinomial <- residuals(nb)
      }, error = function(e) {
        tp$NegativeBinomial <- NULL
      })
    }

    if ("Hurdle" %in% models) {
      tryCatch({
        tp$Hurdle <- residuals(hurdle_mod)
      }, error = function(e) {
        tp$Hurdle <- NULL
      })
    }

    if (any(sapply(tp, is.null))) {
      warning("Error occurred during model fitting. Skipping QQ plot.")
    } else {
      for (i in seq_along(tp)) {
        temp <- data.frame(sam = tp[[i]], Model = models_labels[i])
        df <- rbind(df, temp)
      }

      num_models <- sum(models %in% models_names)
      num_cols <- ifelse(num_models > 2, 2, num_models)
      num_rows <- ceiling(num_models / num_cols)

      gg_list <- list()
      for (i in seq_along(tp)) {
        gg <- df %>%
          filter(Model == models_labels[i]) %>%
          ggplot(aes(sample = sam)) +
          geom_qq(color = meta_colors$model[i]) +
          geom_qq_line(color = 'black') +
          ylab('Residuals') +
          xlab('Theoretical Quantile') +
          ggtitle(paste0(x, ' Q-Q Plot - ', models_labels[i]))

        if (i == num_models) {
          gg <- gg + scale_color_manual(values = meta_colors$model)
        } else {
          gg <- gg + theme(legend.position = 'none')
        }

        gg_list[[i]] <- gg
      }

      plo <- do.call(grid.arrange, c(gg_list, ncol = num_cols))
      setwd(path)
      ggsave(paste0(x, '_qqplot_', Sys.Date(), '.png'), plo, width = 8, height = 6, units = "in")

      plo
    }
  }

  ## Estimate parameter

  if (estimate == T){

    if ("GLM" %in% models) {
      glmout = NULL
      tryCatch({
        glmout = tibble(
          Gene = x,
          Estimate =  as.numeric(summary(mod)$coefficients[2, 1]),
          Std.Error= as.numeric(summary(mod)$coefficients[2, 2]),
          Z_Value=as.numeric(summary(mod)$coefficients[2, 3]),
          P_Value = as.numeric(parameters::p_value(mod)[2,2]),
          CI_lower= Estimate - as.numeric(summary(mod)$coefficients[2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(mod)$coefficients[2, 2])*1.96,
          chisq = lrtest(mod,mod1)[2,4],
          lrt_pval = lrtest(mod,mod1)[2,5],
          AIC = AIC(mod),
          BIC = BIC(mod),
          model = 'glm'
        )
      }, error = function(e) {
        glmout = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          chisq = 'error',
          lrt_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'glm'
        )
      })
    }


    if ("LME" %in% models) {
      lmeout = NULL
      tryCatch({
        lmeout = tibble(
          Gene = x,
          Estimate =  as.numeric(summary(lme_mod)$coefficients[2, 1]),
          Std.Error= as.numeric(summary(lme_mod)$coefficients[2, 2]),
          Z_Value=as.numeric(summary(lme_mod)$coefficients[2, 3]),
          P_Value = as.numeric(parameters::p_value(lme_mod)[2,2]),
          CI_lower= Estimate - as.numeric(summary(lme_mod)$coefficients[2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(lme_mod)$coefficients[2, 2])*1.96,
          chisq = lrtest(lme_mod,lme_mod1)[2,4],
          lrt_pval = lrtest(lme_mod,lme_mod1)[2,5],
          AIC = AIC(lme_mod),
          BIC = BIC(lme_mod),
          model = 'lme'
        )
      }, error = function(e) {
        lmeout = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          chisq = 'error',
          lrt_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'lme'
        )
      })
    }

    if ("ZeroInflated" %in% models) {
      zeroout <- NULL  # Define zeroout outside tryCatch block
      tryCatch({
        zeroout = tibble(
          Gene = x,
          Estimate =  as.numeric(summary(zero_mod)$coefficients$count[2, 1]),
          Std.Error= as.numeric(summary(zero_mod)$coefficients$count[2, 2]),
          Z_Value=as.numeric(summary(zero_mod)$coefficients$count[2, 3]),
          P_Value = as.numeric(parameters::p_value(zero_mod)[2,2]),
          CI_lower= Estimate - as.numeric(summary(zero_mod)$coefficients$count[2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(zero_mod)$coefficients$count[2, 2])*1.96,
          chisq = lrtest(zero_mod,zero_mod1)[2,4],
          lrt_pval = lrtest(zero_mod,zero_mod1)[2,5],
          AIC = AIC(zero_mod),
          BIC= BIC(zero_mod),
          model = 'zero inflated'
        )
      }, error = function(e) {
        zeroout = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          chisq = 'error',
          lrt_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'zero inflated'
        )
      })
    }

    if ("ZeroInflNegBinomial" %in% models) {
      zero_neg_out = NULL
      tryCatch({
        zero_neg_out= tibble(
          Gene = x,
          Estimate =  as.numeric(summary(zero_neg)$coefficients$count[2, 1]),
          Std.Error= as.numeric(summary(zero_neg)$coefficients$count[2, 2]),
          Z_Value=as.numeric(summary(zero_neg)$coefficients$count[2, 3]),
          P_Value = as.numeric(parameters::p_value(zero_neg)[2,2]),
          CI_lower= Estimate - as.numeric(summary(zero_neg)$coefficients[count][2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(zero_neg)$coefficients[count][2, 2])*1.96,
          chisq = lrtest(zero_neg,zero_neg1)[2,4],
          lrt_pval = lrtest(zero_neg,zero_neg1)[2,5],
          AIC = AIC(zero_neg),
          BIC = BIC(zero_neg),
          model = 'zero neg'
        )
      }, error = function(e) {
        zero_neg_out = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          chisq = 'error',
          lrt_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'zero neg'
        )
      })
    }

    if ("NegativeBinomial" %in% models) {
      neg_bin_out = NULL
      tryCatch({
        neg_bin_out = tibble(
          Gene = x,
          Estimate =  as.numeric(summary(nb)$coefficients[2, 1]),
          Std.Error= as.numeric(summary(nb)$coefficients[2, 2]),
          Z_Value=as.numeric(summary(nb)$coefficients[2, 3]),
          P_Value = as.numeric(parameters::p_value(nb)[2,2]),
          CI_lower= Estimate - as.numeric(summary(nb)$coefficients[2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(nb)$coefficients[2, 2])*1.96,
          chisq = lrtest(nb,nb1)[2,4],
          lrt_pval = lrtest(nb,nb1)[2,5],
          AIC = AIC(nb),
          BIC = BIC(nb),
          model = 'nb'
        )
      }, error = function(e) {
        neg_bin_out = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          chisq = 'error',
          lrt_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'nb'
        )
      })
    }

    if ("Hurdle" %in% models) {
      hurdle_mod_out = NULL
      tryCatch({
        hurdle_mod_out = tibble(
          Gene = x,
          Estimate =  as.numeric(summary(hurdle_mod)$coefficients$count[2, 1]),
          Std.Error= as.numeric(summary(hurdle_mod)$coefficients$count[2, 2]),
          Z_Value=as.numeric(summary(hurdle_mod)$coefficients$count[2, 3]),
          P_Value = as.numeric(parameters::p_value(hurdle_mod)[2,2]),
          CI_lower= Estimate - as.numeric(summary(hurdle_mod)$coefficients$count[2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(hurdle_mod)$coefficients$count[2, 2])*1.96,
          chisq = lrtest(hurdle_mod,hurdle_mod1)[2,4],
          lrt_pval = lrtest(hurdle_mod,hurdle_mod1)[2,5],
          AIC = AIC(hurdle_mod),
          BIC= BIC(hurdle_mod),
          model = 'hurdle'
        )
      }, error = function(e) {
        hurdle_mod_out = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          chisq = 'error',
          lrt_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'hurdle'
        )
      })
    }

    if ("Quasipoisson" %in% models) {
      qua_mod_out = NULL
      tryCatch({
        qua_mod_out = tibble(
          Gene = x,
          Estimate =  as.numeric(summary(qua)$coefficients[2, 1]),
          Std.Error= as.numeric(summary(qua)$coefficients[2, 2]),
          Z_Value=as.numeric(summary(qua)$coefficients[2, 3]),
          P_Value = as.numeric(parameters::p_value(qua)[2,2]),
          CI_lower= Estimate - as.numeric(summary(qua)$coefficients[2, 2])*1.96,
          CI_upper = Estimate + as.numeric(summary(qua)$coefficients[2, 2])*1.96,
          deviance = anova(qua,qua1, test= 'F')[2,5],
          F_pval = anova(qua,qua1, test= 'F')[2,6],
          AIC = AIC(qua),
          BIC = BIC(qua),
          model = 'quasi'
        )
      }, error = function(e) {
        qua_mod_out = tibble(
          Gene = x,
          Estimate =  'error',
          Std.Error= 'error',
          Z_Value='error',
          P_Value = 'error',
          CI_lower= 'error',
          CI_upper = 'error',
          deviance = 'error',
          F_pval = 'error',
          AIC = 'error',
          BIC = 'error',
          model = 'quasi'
        )
      })
    }


    return(bind_rows(
      if ("GLM" %in% models) glmout,
      if ("LME" %in% models) lmeout,
      if ("ZeroInflated" %in% models) zeroout,
      if ("ZeroInflNegBinomial" %in% models) zero_neg_out,
      if ("NegativeBinomial" %in% models) neg_bin_out,
      if ("Hurdle" %in% models) hurdle_mod_out,
      if ("Quasipoisson" %in% models) qua_mod_out
    ))
  }



  ## AVG effect plot


  if (log_plot == TRUE) {
    glmout <- NULL
    lmeout <- NULL
    zeroout <- NULL
    zero_neg_out <- NULL
    neg_bin_out <- NULL
    hurdle_mod_out <- NULL
    qua_mod_out <- NULL
    gene  = tibble(Gene = x)

    if ("GLM" %in% models) {
      glmout = NULL

      tryCatch({glmout <- tibble(glm_pred = predict(mod, type = 'response'))
      }, error = function(e){
        glmout <- tibble(glm_pred = 'error')
      })
    }

    if ("LME" %in% models) {
      lmeout <- NULL

      tryCatch({
        lmeout <- tibble(lme_pred = predict(lme_mod, type = 'response'))
      }, error = function(e) {
        lmeout <- tibble(lme_pred = 'error')
      })
    }

    if ("ZeroInflated" %in% models) {
      zeroout <- NULL

      tryCatch({
        zeroout <- tibble(zero_mod_pred = predict(zero_mod, type = 'response'))
      }, error = function(e) {
        zeroout <- tibble(zero_mod_pred = 'error')
      })
    }

    if ("ZeroInflNegBinomial" %in% models) {
      zero_neg_out <- NULL

      tryCatch({
        zero_neg_out <- tibble(zero_neg_pred = predict(zero_neg, type = 'response'))
      }, error = function(e) {
        zero_neg_out <- tibble(zero_neg_pred = 'error')
      })
    }

    if ("NegativeBinomial" %in% models) {
      neg_bin_out <- NULL

      tryCatch({
        neg_bin_out <- tibble(nb_pred = predict(nb, type = 'response'))
      }, error = function(e) {
        neg_bin_out <- tibble(nb_pred = 'error')
      })
    }

    if ("Hurdle" %in% models) {
      hurdle_mod_out <- NULL

      tryCatch({
        hurdle_mod_out <- tibble(hurdle_pred = predict(hurdle_mod, type = 'response'))
      }, error = function(e) {
        hurdle_mod_out <- tibble(hurdle_pred = 'error')
      })
    }

    if ("Quasipoisson" %in% models) {
      qua_mod_out <- tibble(qua_pred = predict(qua, type = 'response'))
    }

    return(bind_cols(
      gene,
      glmout,
      lmeout,
      zeroout,
      zero_neg_out,
      neg_bin_out,
      hurdle_mod_out,
      qua_mod_out
    ))
  }


  ## Rootogram

  if (rootogram == TRUE) {
    plotList <- list()

    if ("GLM" %in% models) {
      a <- rootogram(mod, plot = FALSE)
      p1 <- autoplot(a) + labs(title = paste(x, ' - GLM'))
      plotList <- c(plotList, list(p1))
    }

    if ("ZeroInflNegBinomial" %in% models) {
      p <- rootogram(zero_neg, plot = FALSE, sub = 'Zero Negative')
      p2 <- autoplot(p) + labs(title = paste(x, ' - Zero Inflated Negative Binomial'))
      plotList <- c(plotList, list(p2))
    }

    if ("NegativeBinomial" %in% models) {
      c <- rootogram(nb, plot = FALSE)
      p3 <- autoplot(c) + labs(title = paste(x, ' - Negative Binomial'))
      plotList <- c(plotList, list(p3))
    }

    if ("ZeroInflated" %in% models) {
      z <- rootogram(zero_mod, plot = FALSE)
      p4 <- autoplot(z) + labs(title = paste(x, ' - Zero Inflated'))
      plotList <- c(plotList, list(p4))
    }

    if ("Hurdle" %in% models) {
      h <- rootogram(hurdle_mod, plot = FALSE, width = 1)
      p5 <- autoplot(h) + labs(title = paste(x, ' - Hurdle'))
      plotList <- c(plotList, list(p5))
    }

    if (length(plotList) > 0) {
      plo <- do.call(grid.arrange, c(plotList, ncol = 1))

      setwd(path)
      ggsave(paste0(x, '_rootogram_', Sys.Date(), '.pdf'), plo, units = 'in', width = 3, height = 6)
      plo
    }
  }
}
