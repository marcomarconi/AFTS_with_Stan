# some helper functions we'll use throughout

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(loglik, ids = NULL) {
  if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
  rowSums(loglik)
}

# for printing comparisons later
rbind_print <- function(...) {
  round(rbind(...), digits = 2)
}


plot_ks <- function(ks, ids, thres = 0.6) {
  dat_ks <- data.frame(ks = ks, ids = ids)
  ggplot(dat_ks, aes(x = ids, y = ks)) +
    geom_point(aes(color = ks > thres), shape = 3, show.legend = FALSE) +
    geom_hline(yintercept = thres, linetype = 2, color = "red2") +
    scale_color_manual(values = c("cornflowerblue", "darkblue")) +
    labs(x = "Data point", y = "Pareto k") +
    ylim(-0.5, 1.5)
}


k_thres <- 0.7




LOO_LFO <- function(fit, df, L, M) {
  N <- nrow(df)
  approx_elpds_Msap <- rep(NA, N)
  
  # initialize the process for i = L
  past <- 1:L
  oos <- (L + 1):(L + M)
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]
  fit_past <- update(fit, newdata = df_past, recompile = FALSE, cores=4)
  loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
  loglikm <- rowSums(loglik[, oos, drop=F])
  approx_elpds_1sap[L + 1] <- log_mean_exp(loglikm)
  
  # iterate over i > L
  i_refit <- L
  refits <- L
  ks <- NULL
  for (i in (L + 1):(N - M)) {
    past <- 1:i
    oos <- (i + 1):(i + M)
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
    
    logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
    psis_obj <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_obj)
    ks <- c(ks, k)
    if (k > k_thres) {
      # refit the model based on the first i observations
      i_refit <- i
      refits <- c(refits, i)
      fit_past <- update(fit_past, newdata = df_past, recompile = FALSE)
      loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
      loglikm <- rowSums(loglik[, oos, drop=F])
      approx_elpds_Msap[i + 1] <- log_mean_exp(loglikm)
    } else {
      lw <- weights(psis_obj, normalize = TRUE)[, 1]
      loglikm <- rowSums(loglik[, oos, drop=F])
      approx_elpds_Msap[i + 1] <- log_sum_exp(lw + loglikm)
    }
  }
  
  return(c("ks"=ks, "approx_elpds_Msap"=approx_elpds_Msap))
}


LOO_LFO_1 <- function(fit, y, L, model_data=list() ,SEED=1983) {
  N = length(y)
  approx_elpds_1sap <- rep(NA, N)

  # initialize the process for i = L
  past <- 1:L
  oos <- L + 1
  y_past <- y[past]
  y_oos <- y[oos:length(y)]
  fit_past <- stan(fit=fit,  data=c(list(N = length(y_past), y = y_past, M = length(y)-L, m =  y_oos), model_data) , chains=2, cores=2,  iter=1000, seed = SEED)
  loglik <- extract_log_lik(fit_past)
  approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])
  
  # iterate over i > L
  i_refit <- L
  refits <- L
  ks <- NULL
  for (i in (L + 1):(N - 1)) {
    print(i)
    past <- 1:i
    oos <- i + 1
    y_past <- y[past]
    y_oos <- y[oos:length(y)]
    logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
    psis_obj <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_obj)
    ks <- c(ks, k)
    if (k > k_thres) {
      # refit the model based on the first i observations
      i_refit <- i
      refits <- c(refits, i)
      fit_past <- stan(fit=fit,  data=c(list(N = length(y_past), y = y_past, K = 1, M = length(y_oos), m =  y_oos), model_data), chains=2, cores=2,  iter=1000, seed = SEED)
      loglik <- extract_log_lik(fit_past)
      approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])
    } else {
      lw <- weights(psis_obj, normalize = TRUE)[, 1]
      approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
    }
  }
  return(list(approx_elpds_1sap=approx_elpds_1sap, ks=ks, refits=refits))
}

