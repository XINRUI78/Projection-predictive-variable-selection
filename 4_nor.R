library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)

# Run projpred for one prior
run_nor <- function(i, ndev, nval, n.para, beta0, beta,
                          nterms_max = 30,
                          ns = 2000) {
  set.seed(i)
  SEED <- as.integer(i)
  
  data.dev <- generate_ss(ndev, n.para, beta0, beta)
  data.val <- generate_ss(nval, n.para, beta0, beta)
  
  xval <- data.val[, -1]
  yval <- data.val[, 1]
  
  # Choose the full model as the reference model
  all_vars <- colnames(data.dev)[-1]
  ref_formula <- as.formula(
    paste("y", paste(all_vars, collapse = " + "), sep = " ~ ")
  )
  
  # Fit Bayesian reference model
  ref_fit <- stan_glm(
    formula = ref_formula,
    data = data.dev,
    family = binomial(link = "logit"),
    prior = normal(0, 1),
    prior_intercept = normal(0, 1),
    QR = TRUE,
    seed = SEED,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    chains = 4
  )
  ref_fit$call$formula <- ref_formula
  ref_fit$call$data <- data.dev
  ref_fit$call$prior <- prior_obj
  ref_fit$call$seed <- SEED
  
  # Reference model predictions
  p_ref <- colMeans(
    posterior_epred(ref_fit, newdata = xval)
  )
  
  ref_row <- c(
    0.3, 0.8, ndev, "ref-nor",
    measures(yval, p_ref),
    rep(1, n.para),
    NA
  )
  
  # Projection predictive variable selection
  vs <- cv_varsel(
    ref_fit,
    method = "forward",
    cv_method = "LOO",
    validate_search = TRUE,
    nterms_max = min(nterms_max, n.para),
    seed = SEED + 1000
  )
  
  # Obtain a vector that indicates whether each variable is selected or not, with the suggested size
  get_selected <- function(vs, k) {
    sel_vars <- ranking(vs)$fulldata[1:k]
    as.integer(all_vars %in% sel_vars)
  }
  size_suggested <- suggest_size(vs)
  
  # Obtain a vector that indicates whether each variable is selected or not, with the best ELPD size
  get_best_size <- function(vs) {
    perf <- performances(vs)$submodels
    perf$size[which.max(replace(perf$elpd, is.na(perf$elpd), -Inf))]
  }
  size_best <- get_best_size(vs)
  
  # Project to suggested size
  proj_suggest <- project(
    vs,
    nterms = size_suggested,
    seed = SEED + 10000,
    ns = ns
  )
  
  # Project to best ELPD size
  proj_best <- project(
    vs,
    nterms = size_best,
    seed = SEED + 10001,
    ns = ns
  )
  
  # Predictions from projected models
  p_suggest <- as.vector(
    proj_linpred(
      proj_suggest,
      newdata = xval,
      integrated = TRUE,
      transform = TRUE
    )$pred
  )
  
  p_best <- as.vector(
    proj_linpred(
      proj_best,
      newdata = xval,
      integrated = TRUE,
      transform = TRUE
    )$pred
  )
  
  suggest_row <- c(
    0.3, 0.8, ndev, "nor-1se",
    measures(yval, p_suggest),
    get_selected(vs, size_suggested),
    NA
  )
  
  best_row <- c(
    0.3, 0.8, ndev, "nor-best",
    measures(yval, p_best),
    get_selected(vs, size_best),
    NA
  )
  
  rbind(ref_row, suggest_row, best_row)
}

