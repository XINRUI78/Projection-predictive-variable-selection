library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)

# Run projpred for one prior
run_one_prior <- function(data.dev, data.val, ndev, n.para,
                          prior_obj, prior_name,
                          seed_fit, seed_cv, seed_proj,
                          nterms_max = 30,
                          K = 5,
                          ns = 2000) {

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
    prior = prior_obj,
    prior_intercept = normal(0, 1),
    QR = TRUE,
    seed = seed_fit,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    chains = 4
  )
  ref_fit$call$formula <- ref_formula
  ref_fit$call$data <- data.dev
  ref_fit$call$prior <- prior_obj
  ref_fit$call$seed <- seed_fit

  # Reference model predictions
  p_ref <- colMeans(
    posterior_epred(ref_fit, newdata = xval)
  )

  ref_row <- c(
    0.3, 0.8, ndev, paste0("ref-", prior_name),
    measures(yval, p_ref),
    rep(1, n.para),
    NA
  )

  # Projection predictive variable selection
  vs <- cv_varsel(
      ref_fit,
      method = "forward",
      cv_method = "kfold",
      validate_search = TRUE,
      K = K,
      nterms_max = min(nterms_max, n.para),
      parallel = TRUE,
      seed = seed_cv
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
    seed = seed_proj,
    ns = ns
  )

  # Project to best ELPD size
  proj_best <- project(
    vs,
    nterms = size_best,
    seed = seed_proj + 1,
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
    0.3, 0.8, ndev, paste0(prior_name, "-1se"),
    measures(yval, p_suggest),
    get_selected(vs, size_suggested),
    NA
  )

  best_row <- c(
    0.3, 0.8, ndev, paste0(prior_name, "-best"),
    measures(yval, p_best),
    get_selected(vs, size_best),
    NA
  )

  rbind(ref_row, suggest_row, best_row)
}

# Run one simulation replicate
run_projpred_one <- function(i, ndev, n.para, beta0, beta, nval) {

  set.seed(i)
  SEED <- as.integer(i)

  data.dev <- generate_ss(ndev, n.para, beta0, beta)
  data.val <- generate_ss(nval, n.para, beta0, beta)

  res_norm <- run_one_prior(
    data.dev = data.dev,
    data.val = data.val,
    ndev = ndev,
    n.para = n.para,
    prior_obj = normal(0, 1),
    prior_name = "nor",
    seed_fit = SEED,
    seed_cv = SEED + 1000,
    seed_proj = SEED + 10000,
    nterms_max = 25,
    K = 5,
    ns = 2000
  )

  res_laplace <- run_one_prior(
    data.dev = data.dev,
    data.val = data.val,
    ndev = ndev,
    n.para = n.para,
    prior_obj = laplace(location = 0, scale = 1),
    prior_name = "la",
    seed_fit = SEED + 2000,
    seed_cv = SEED + 3000,
    seed_proj = SEED + 40000,
    nterms_max = 25,
    K = 5,
    ns = 2000
  )

  method_result <- rbind(res_norm, res_laplace)

  colnames(method_result) <- c(
    "prevalence",
    "anticipated c-stat",
    "ndev",
    "method",
    "calibration slope",
    "calibration in the large",
    "auc",
    "rmspe",
    paste0("varsel", 1:n.para),
    "option"
  )

  method_result
}
