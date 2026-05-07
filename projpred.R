set.seed(123)
## Setup
#### Data
# development dataset of size 752
data.dev <- read.csv("data_dev.csv")
sam <- sample(1:752, 376) 
data.dev <- data.dev[sam,]
x <- data.dev[,-1]
y <- data.dev[,1]

# validation dataset of size 10000
data.val <- read.csv("data_val.csv")
xval <- data.val[,-1]
yval <- data.val[,1]

#### Packages
library(tidyverse)
library(rstanarm)
library(brms)
library(projpred)
library(tidybayes)
library(broom.mixed)
library(pROC)
library(patchwork)
library(cutpointr)
library(gt)
library(coda)
library(loo)
library(glmnet)
library(detectseparation)
library(brglm2)
library(foreach)
library(MASS)
options(mc.cores = parallel::detectCores())

SEED = 7
set.seed(SEED)
refform <- as.formula(paste("y", paste(colnames(data.dev)[-1], collapse=" + "), sep=" ~ "))
ref_fit <- stan_glm(formula = refform, 
                    data = data.dev,
                    family = binomial(link = "logit"), 
                    prior = laplace(location = 0, scale = 1), 
                    prior_intercept = normal(0, 1), 
                    QR = TRUE, 
                    seed = SEED, 
                    adapt_delta = 0.99,
                    iter = 4000,
                    cores = 4,
                    chains = 4)

## Check model diagnostics

summary(ref_fit) 
pp_check(ref_fit, plotfun = "bars", nreps = 500)
loo_reference_fit <- loo(ref_fit)
loo_reference_fit
plot(loo_reference_fit)

## Projection predictive variable selection
# full data search
varsel_pro <- cv_varsel(ref_fit, 
                        method = 'forward', 
                        validate_search = F,
                        cores = 8,
                        seed = 11)
(vif <- ranking(varsel_pro)$fulldata)
(size_suggested <- suggest_size(varsel_pro))

library(pander)
s <- summary(
  varsel_pro,
  stats = c("elpd", "acc", "auc", "rmse"),
  deltas = FALSE
)
performances(s)$submodels %>%
  dplyr::select(size, elpd, elpd.se, elpd.diff, elpd.diff.se, acc, acc.se, auc, auc.se, rmse, rmse.se) %>% 
  pander::pander(
    split.cell = 80,
    split.table = Inf,
    justify = "left",
    caption = "Feature selection trajectory statistics for the projection submodels"
  )

perf <- performances(varsel_pro)$submodels
size_best <- perf$size[which.max(perf$elpd)]
plot(varsel_pro, stats = c("elpd", "acc", "auc", "rmse")) +
  geom_vline(xintercept = size_suggested, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = size_best, color = "red", linetype = "dashed") 

# 10-fold cross-validation search
varsel_proj <- cv_varsel(ref_fit, 
                         method = 'forward', 
                         cv_method = "kfold",
                         validate_search = T,
                         K = 10,
                         cores = 8,
                         seed = 11)

(vif <- ranking(varsel_proj)$fulldata)
(size_suggested <- suggest_size(varsel_proj))

library(pander)
s <- summary(
  varsel_proj,
  stats = c("elpd", "acc", "auc", "rmse"),
  deltas = FALSE
)
performances(s)$submodels %>%
  dplyr::select(size, elpd, elpd.se, elpd.diff, elpd.diff.se, acc, acc.se, auc, auc.se, rmse, rmse.se) %>% 
  pander::pander(
    split.cell = 80,
    split.table = Inf,
    justify = "left",
    caption = "Feature selection trajectory statistics for the projection submodels"
  )

perf <- performances(varsel_proj)$submodels
size_best <- perf$size[which.max(perf$elpd)]
plot(varsel_proj, stats = c("elpd", "acc", "auc", "rmse")) +
  geom_vline(xintercept = size_suggested, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = size_best, color = "red", linetype = "dashed") 

## Projection

proj_suggest <- project(varsel_proj, 
                        nterms = size_suggested, 
                        seed = SEED, ns = 2000)
proj_best <- project(varsel_proj, 
                     nterms = size_best, 
                     seed = SEED, ns = 2000)

pred_suggest <- proj_linpred(proj_suggest,
                             newdata = data.val[,-1], 
                             integrated = TRUE,
                             transform = TRUE
)
pred_best <- proj_linpred(proj_best,
                          newdata = data.val[,-1], 
                          integrated = TRUE,
                          transform = TRUE
)

## Performance

measures <- function(yval, p_val) {
  
  eta_val <- log(p_val/(1 - p_val))
  
  # Calibration slope
  fitcal <- speedglm::speedglm(yval ~ eta_val, family = binomial())
  cal_slope <- as.vector(coef(fitcal)[2]) 
  
  # Calibration in the large
  off <- speedglm::speedglm(yval ~ 1, offset = eta_val, family = binomial())
  cal_large <- as.vector(coef(off))
  
  # AUC
  cstat <- pROC::roc(response = yval, predictor = as.vector(p_val), levels = c(0, 1), direction = "<")
  auc <- as.vector(cstat$auc)
  
  # Root mean square prediction error (RMSPE)
  rmspe <- sqrt(mean((p_val - yval)^2))
  return(c(cal_slope, cal_large, auc, rmspe))
}

p_suggest <- pred_suggest$pred[1,]
meas_suggest <- measures(yval, p_suggest)
suggest_stat <- c(size_suggested, meas_suggest)
cat("calibration slope, calibration in the large, auc, rmse values:", meas_suggest)

p_best <- pred_best$pred[1,]
meas_best <- measures(yval, p_best)
best_stat <- c(size_best, meas_best)
cat("calibration slope, calibration in the large, auc, rmse values:", meas_best)

backward_pvalue <- function(data, response, p_threshold = 0.05) {
  
  # Start with the full model for backward elimination
  current_formula <- as.formula(paste(response, "~ ."))
  model <- glm(current_formula, data = data, family = binomial)
  predictors <- setdiff(names(data), response)
  
  # Main loop for stepwise selection
  while (TRUE) {
    removed <- FALSE
    
    # Backward Step: Remove predictors with p-value > threshold
    p_values <- summary(model)$coefficients[-1, 4]  # Exclude intercept
    names(p_values) <- rownames(summary(model)$coefficients)[-1]  # Get predictor names
    max_pval <- max(p_values, na.rm = TRUE)
    
    # If there is a predictor to remove with p-value > threshold
    if (max_pval > p_threshold) {
      worst_predictor <- names(which.max(p_values))
      current_formula <- as.formula(paste(response, "~", paste(setdiff(predictors, worst_predictor), collapse = " + ")))
      model <- glm(current_formula, data = data, family = binomial)
      predictors <- setdiff(predictors, worst_predictor)
      removed <- TRUE
    }
    
    
    # Stop if no further changes can be made
    if (!removed) {
      break
    }
  }
  
  return(model)
}
backward_lrt <- function(data, response, p_threshold = 0.05) {
  predictors <- setdiff(names(data), response)
  current_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  current_model <- glm(current_formula, data = data, family = binomial)
  
  while (length(predictors) > 1) {
    p_values <- sapply(predictors, function(var) {
      reduced_predictors <- setdiff(predictors, var)
      reduced_formula <- as.formula(paste(response, "~", paste(reduced_predictors, collapse = " + ")))
      reduced_model <- glm(reduced_formula, data = data, family = binomial)
      
      # Likelihood Ratio Test
      lrt <- anova(reduced_model, current_model, test = "LRT")
      pval <- lrt[2, "Pr(>Chi)"]
      return(pval)
    })
    
    max_p <- max(p_values, na.rm = TRUE)
    if (max_p > p_threshold) {
      worst_var <- names(which.max(p_values))
      predictors <- setdiff(predictors, worst_var)
      current_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
      current_model <- glm(current_formula, data = data, family = binomial)
    } else {
      break
    }
  }
  
  return(current_model)
}
# Main function
back_logit <- function(data, response, crit, p_threshold = 0.05) {
  data[[response]] <- as.factor(data[[response]])
  
  full_model <- glm(as.formula(paste(response, "~ .")), data = data, family = binomial)
  
  if (crit == "AIC" || crit == "BIC") {
    final_model <- step(full_model, direction = "backward", k = ifelse(crit == "AIC", 2, log(nrow(data))))
    
  } else if (crit == "p_value") {
    final_model <- backward_pvalue(data, response, p_threshold)
    
  } else if (crit == "LRT") {
    final_model <- backward_lrt(data, response, p_threshold)
    
  } else {
    stop("Invalid criterion. Use 'AIC', 'BIC', 'p_value', or 'LRT'.")
  }
  
  back_predictors <- names(coef(final_model))[-1]  # Exclude intercept
  all_predictors <- setdiff(names(data), response)
  indicator <- ifelse(all_predictors %in% back_predictors, 1, 0)
  
  return(list(
    varsel_back = indicator,
    backmodel = final_model
  ))
}


unilogit <- function(x, y, p_threshold = 0.05) {
  
  # Ensure y is a factor for logistic regression
  #y <- as.factor(y)
  
  # Initialize storage for results
  results <- list()
  varsel_uni <- numeric(ncol(x))  # Binary row for significance
  
  # Loop through each column by number
  for (i in 1:ncol(x)) {
    
    # Fit univariate logistic regression using column index
    predictor <- x[, i]
    model <- glm(y ~ predictor, family = binomial)
    
    # Extract p-value
    summary_stats <- summary(model)$coefficients
    p_value <- summary_stats[2, 4]
    
    # Determine significance
    varsel_uni[i] <- ifelse(p_value < p_threshold, 1, 0)
  }  
  
  # Create the final model using significant variables
  significant_vars <- which(varsel_uni == 1)
  if (length(significant_vars) > 0) {
    final_formula <- as.formula(paste("y ~", paste(paste0("x[,", significant_vars, "]"), collapse = " + ")))
    final_model <- glm(final_formula, family = binomial)
  } else {
    final_model <- glm(y ~ 1, family = binomial)
  }
  
  # Return results
  return(list(
    varsel_uni = varsel_uni,
    uni = final_model
  ))
}


# Inputs

#' @param x        a matrix, the design matrix
#' @param y        a vector, the binary response variable
#' @param bn       a scalar, the number of bootstrap iterations
#' @param method  a string, values "ridge" or "lasso"Numeric where to trim (proportion)
#' @param nfolds   a scalar, the number of cross-validation folds
#' @param parallel a logical, TRUE for parallel computing
#' @param boot     a logical, TRUE means sample with replacement (the default)
#'
#' Outputs: a list with two elements.
#'          - beta.boot is the calculated vector of regression coefficients
#'          - lambda.boot is vector of beta's


mod_penal_ave_foreach <- function(x, y, bn=20, method="ridge", nfolds=10, f = nfolds/(nfolds-1)-1,
                                  parallel=TRUE,  boot=TRUE){
  library(glmnet)
  library(foreach)
  
  if (method=="lasso") a=1
  if (method=="ridge") a=0
  
  pred        <- NULL ; beta.boot   <- NULL
  lambda.boot <- NULL ; cv.boot     <- NULL
  yval        <- y
  xval        <- x
  npred       <- ncol(x)
  data        <- cbind(x,y)
  
  fit <- glmnet(x, y, type.measure = "deviance", family = "binomial", alpha = a, standardize = TRUE, parallel = T)
  
  # Obtain the default sequence of lambda, and include some smaller values too
  lambdaseq <- fit$lambda
  lambdaseq <- unique(c(lambdaseq, seq( lambdaseq[length(lambdaseq)], lambdaseq[length(lambdaseq)]/20,
                                        -lambdaseq[length(lambdaseq)]/20 )))
  
  # cvd is the cross-validated deviance obtained in each of the 'over-sampled datasets'
  cvd <- NULL
  
  cv.boot <- foreach(i=1:bn, .combine='rbind', .packages=c('MASS','RcppNumerical', 'pROC', 'brglm2', 'detectseparation',
                                                           'RcppNumerical', 'logistf', 'glmnet')) %do% {
                                                             
                                                             
                                                             withwarnings <- function(expr) {
                                                               
                                                               val        <- NULL
                                                               myWarnings <- NULL
                                                               wHandler   <- function(w) {
                                                                 myWarnings <<- c(myWarnings, w$message)
                                                                 invokeRestart("muffleWarning")
                                                               }
                                                               
                                                               myError  <- NULL
                                                               eHandler <- function(e) {
                                                                 myError <<- e$message
                                                                 NULL
                                                               }
                                                               val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
                                                               list(value = val, warnings = myWarnings, error=myError)
                                                             }
                                                             
                                                             bs     <- sample(nrow(data), replace=T)
                                                             databs <- data[bs,]
                                                             
                                                             if (boot==FALSE) databs<-data
                                                             
                                                             if (f==0) databs  <- data else
                                                             {
                                                               
                                                               bs2     <- sample(nrow(data)*f, replace=T)
                                                               databs2 <- data[bs2,]
                                                               databs  <- rbind(databs,databs2)
                                                             }
                                                             
                                                             xbs   <- databs[,1:npred]
                                                             ybs   <- databs[,npred+1]
                                                             
                                                             fitbs <- withwarnings(cv.glmnet(xbs, ybs, type.measure = "deviance", family = "binomial", alpha = a, nfolds = nfolds,
                                                                                             standardize = TRUE, parallel = parallel, lambda = lambdaseq))
                                                             
                                                             
                                                             if ( length(fitbs$error) == 0 ) {
                                                               fitbs  <- fitbs$value
                                                               cvd    <- fitbs$cvm
                                                               
                                                             } else
                                                               
                                                             {
                                                               cvd <- cvd
                                                             }
                                                             
                                                             rbind(cvd)
                                                             
                                                           }
  
  b <- cbind(lambdaseq, colMeans(cv.boot, na.rm = TRUE))
  
  # Selected the value of lambda that minimises the cross-validated deviance
  
  lambda.boot <- b[order(b[,2]),] [1,1]
  
  
  # Obtain the solution using ridge/lasso at the selected value of lambda for the *original data*
  
  fit <- glmnet(x, y, type.measure = "deviance", family = "binomial", alpha = a, standardize = TRUE,
                parallel = T, lambda = lambdaseq)
  
  beta.boot   <- as.vector(coef(fit, s = lambda.boot))
  
  return(list("beta.boot" = beta.boot, "lambda.boot" = lambda.boot))
}

# frequentist full model
fullfit <- glm(y ~ ., data = data.dev, family = 'binomial')
eta_full <- as.matrix(cbind(1,xval))%*%coef(fullfit)
p_full <- as.vector(1/(1+exp(-eta_full)))
full_stat <- c(30, measures(yval, p_full))

# BE with p-value threshold of 0.05
p_threshold = 0.05
back_05 <- back_logit(data.dev, "y", "LRT", 0.05) 
varsel_back_05 <- back_05$varsel_back
backmodel_05 <- back_05$backmodel
back_eta_05 <- as.matrix(cbind(1,xval[,varsel_back_05 == 1]))%*%coef(backmodel_05)
back_p_05 <- as.vector(1/(1+exp(-back_eta_05)))
back_stat <- c(sum(varsel_back_05), measures(yval, back_p_05))

# UVS with p-value threshold of 0.05
p_threshold <- 0.05
uni_05 <- unilogit(x, y, p_threshold) 
varsel_uni_05 <- uni_05$varsel_uni
unimodel_05 <- uni_05$uni
uni_eta_05 <- as.matrix(cbind(1, xval[, varsel_uni_05 == 1])) %*% coef(unimodel_05)
uni_p_05 <- as.vector(1 / (1 + exp(-uni_eta_05)))
uni_stat <- c(sum(varsel_uni_05), measures(yval, uni_p_05))

# LASSO-min
lasso <- glmnet::cv.glmnet(as.matrix(x), y, alpha = 1, family = "binomial", type.measure = "deviance")
lambda_min <- lasso$lambda.min 
lassomin_p <- as.vector(predict(lasso, as.matrix(xval), s = lambda_min, type="response"))
varsel_min <- ifelse(as.numeric(coef(lasso, s = lambda_min)[-1]) != 0, 1, 0)
LSmin_stat <- c(sum(varsel_min), measures(yval, lassomin_p))

# LASSO-mod
nfolds <- 10
f <- nfolds / (nfolds - 1) - 1
mod_lasso <- mod_penal_ave_foreach(x = as.matrix(x), y = y, bn = 20, method = "lasso", f = f,
                                   parallel = FALSE, nfolds = nfolds, boot = TRUE)
eta_mod_lasso <- as.matrix(cbind(1,xval))%*% mod_lasso$beta.boot
p_mod_lasso <- as.vector(1/(1+exp(-eta_mod_lasso)))
varsel_mod_lasso <- ifelse(as.numeric(mod_lasso$beta.boot)[-1] != 0, 1, 0)
LSmod_stat <- c(sum(varsel_mod_lasso), measures(yval, p_mod_lasso))

# Bayesian reference model
p_ref <- colMeans(posterior_epred(ref_fit, newdata = xval))
ref_stat <- c(30, measures(yval, p_ref))

model <- c('Bayesian reference model',
           'Projected submodel (suggest)',
           'Projected submodel (best)',
           'full model',
           'BE-0.05',
           'UVS-0.05',
           'LASSO-min',
           'LASSO-mod')
df <- data.frame(
  Model = model,
  do.call(rbind, list(
    ref_stat,
    suggest_stat,
    best_stat,
    full_stat,
    back_stat,
    uni_stat,
    LSmin_stat,
    LSmod_stat))
)
colnames(df)[-1] <- c("Number of predictors","Calibration slope", "Calibration in the large", "AUC", "RMSE")
df[-1] <- round(df[-1], 3)
pander(
  df,
  justify = c("left", rep("center", ncol(df) - 1)),
  caption = "Model performance comparison"
)
