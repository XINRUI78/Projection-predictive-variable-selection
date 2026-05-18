library(doParallel)
library(foreach)

freq_n <- perform_n(ndev, n.para, beta0, beta, nval)
freq_n_2 <- perform_n(ndev1, n.para, beta0, beta, nval)
freq_n_4 <- perform_n(ndev2, n.para, beta0, beta, nval)

perform_n <- function(ndev, n.para, beta0, beta, nval){

  # Register parallel backend with the number of cores available
  n.cores <- parallel::detectCores() - 2  # use one less than total cores to avoid overloading
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)

  # Generate a large dataset (200,000) using function generate_ss
  # Check prevalence and C-statistic;
  #n <- 200000;
  #data <- generate_ss(n, n.para, beta0, beta)
  prev <- 0.3 #prev <- round(mean(data[,1]),2)
  #X <- as.matrix(data[,-1])
  #eta <- rep(beta0, n) + X%*%beta
  #p <- 1/(1+exp(-eta))
  #cstat <- pROC::roc(response = as.numeric(data[,1]), predictor = as.vector(p), levels = c(0, 1), direction = "<")
  auc <- 0.8 #auc <- as.numeric(round(as.vector(cstat$auc),2))

  # create a matrix for each method
  n.loop <- 100 # first try 5 loops to check the code, then 100
  results <- foreach(i = 1: n.loop, .combine = rbind, .packages = c("glmnet", "pROC", "stats", "tidyverse", "rstanarm", "projpred", "loo", "pander"),
                     .export = c("generate_ss", "measures", "back_logit", "backward_pvalue", "mod_penal_ave_foreach", "unilogit")) %dopar%{

                       set.seed(i)

                       # generate a development dataset of size ndev=2000
                       data.dev <- generate_ss(ndev, n.para, beta0, beta)
                       x <- data.dev[,-1]
                       y <- data.dev[,1]

                       # generate a validation dataset of size nval
                       data.val <- generate_ss(nval, n.para, beta0, beta)
                       xval <- data.val[,-1]
                       yval <- data.val[,1]

                       ##########################################################

                       # Model fitting
                       # Initialize matrices for the different methods for this iteration

                       method_result <- matrix(NA, nrow = 7, ncol = 9 + n.para)
                       method_name <- NA
                       SEED <- as.integer(i)


                       # method = 1 for full model
                       fit <- glm(y ~ ., data = data.dev, family = 'binomial')
                       eta_val <- as.matrix(cbind(1,xval))%*%coef(fit)
                       p_val <- as.vector(1/(1+exp(-eta_val)))
                       method_result[1,] <- c(prev, auc, ndev, "MLE", measures(yval, p_val), rep(1, n.para), NA)

                       tryCatch({
                         method_name <- "univariable logistic (p < 0.05)"
                         # method = 3 for univariable logistic (p < 0.05)
                         p_threshold <- 0.05
                         unisum_05 <- unilogit(x, y, p_threshold)
                         varsel_uni_05 <- unisum_05$varsel_uni
                         unimodel_05 <- unisum_05$uni
                         uni_eta_05 <- as.matrix(cbind(1, xval[, varsel_uni_05 == 1])) %*% coef(unimodel_05)
                         uni_p_05 <- as.vector(1 / (1 + exp(-uni_eta_05)))
                         method_result[2,] <- c(prev, auc, ndev, "UVS-5%", measures(yval, uni_p_05), varsel_uni_05, p_threshold)
                       }, error = function(e) {
                         cat(sprintf("Error in %s: %s\n", method_name, e$message))
                       })

                       tryCatch({
                         method_name <- "univariable logistic (p < 0.15)"
                         # method = 4 for univariable logistic (p < 0.15)
                         p_threshold <- 0.15
                         unisum_15 <- unilogit(x, y, p_threshold)
                         varsel_uni_15 <- unisum_15$varsel_uni
                         unimodel_15 <- unisum_15$uni
                         uni_eta_15 <- as.matrix(cbind(1, xval[, varsel_uni_15 == 1])) %*% coef(unimodel_15)
                         uni_p_15 <- as.vector(1 / (1 + exp(-uni_eta_15)))
                         method_result[3, ] <- c(prev, auc, ndev, "UVS-15%", measures(yval, uni_p_15), varsel_uni_15, p_threshold)
                       }, error = function(e) {
                         cat(sprintf("Error in %s: %s\n", method_name, e$message))
                       })

                       tryCatch({
                         method_name <- "Backward Elimination (p < 0.05)"
                         # method = 1 for backward elimination for p value of 0.05
                         p_threshold = 0.05
                         back_05 <- back_logit(data.dev, "y", "p_value", p_threshold)
                         varsel_back_05 <- back_05$varsel_back
                         backmodel_05 <- back_05$backmodel
                         back_eta_05 <- as.matrix(cbind(1,xval[,varsel_back_05 == 1]))%*%coef(backmodel_05)
                         back_p_05 <- as.vector(1/(1+exp(-back_eta_05)))
                         method_result[4,] <- c(prev, auc, ndev, "BE-5%", measures(yval, back_p_05), varsel_back_05, p_threshold)
                       }, error = function(e) {
                         cat(sprintf("Error in %s: %s\n", method_name, e$message))
                       })

                       tryCatch({
                         method_name <- "Backward Elimination (p < 0.15)"
                         # method = 2 for backward elimination for p value of 0.15
                         p_threshold = 0.15
                         back_15 <- back_logit(data.dev, "y", "p_value", p_threshold)
                         varsel_back_15 <- back_15$varsel_back
                         backmodel_15 <- back_15$backmodel
                         back_eta_15 <- as.matrix(cbind(1,xval[,varsel_back_15 == 1]))%*%coef(backmodel_15)
                         back_p_15 <- as.vector(1/(1+exp(-back_eta_15)))
                         method_result[5,] <- c(prev, auc, ndev, "BE-15%", measures(yval, back_p_15), varsel_back_15, p_threshold)
                       }, error = function(e) {
                         cat(sprintf("Error in %s: %s\n", method_name, e$message))
                       })

                       tryCatch({
                         method_name <- "LASSO"
                         # method = 7 for LASSO using lambda.min, method = 8 for LASSO using lambda.1se
                         lasso <- glmnet::cv.glmnet(as.matrix(x), y, alpha = 1, family = "binomial", type.measure = "deviance")
                         lambda_min <- lasso$lambda.min
                         #lambda_1se <- lasso$lambda.1se

                         # validate the fitted models on the validation dataset
                         lassomin_p <- as.vector(predict(lasso, as.matrix(xval), s = lambda_min, type="response"))
                         #lasso1se_p <- as.vector(predict(lasso, as.matrix(xval), s = lambda_1se, type="response"))

                         varsel_min <- ifelse(as.numeric(coef(lasso, s = lambda_min)[-1]) != 0, 1, 0)
                         #varsel_1se <- ifelse(as.numeric(coef(lasso, s = lambda_1se)[-1]) != 0, 1, 0)

                         method_result[6,] <- c(prev, auc, ndev, "LS-min", measures(yval, lassomin_p), varsel_min, lambda_min)
                         #method_result[7,] <- c(prev, auc, ndev, 8, measures(yval, lasso1se_p), varsel_1se, lambda_1se)
                       }, error = function(e) {
                         cat(sprintf("Error in %s: %s\n", method_name, e$message))
                       })

                       tryCatch({
                         method_name <- "modified LASSO"
                         # method = 13 for modified LASSO
                         nfolds <- 10
                         f <- nfolds / (nfolds - 1) - 1
                         mod_lasso <- mod_penal_ave_foreach(x = as.matrix(x), y = y, bn = 20, method = "lasso", f = f,
                                                            parallel = FALSE, nfolds = nfolds, boot = TRUE)
                         eta_mod_lasso <- as.matrix(cbind(1,xval))%*% mod_lasso$beta.boot
                         p_mod_lasso <- as.vector(1/(1+exp(-eta_mod_lasso)))
                         varsel_mod_lasso <- ifelse(as.numeric(mod_lasso$beta.boot)[-1] != 0, 1, 0)
                         method_result[7, ] <- c(prev, auc, ndev, "LS-mod", measures(yval, p_mod_lasso), varsel_mod_lasso, as.numeric(mod_lasso$lambda.boot))
                       }, error = function(e) {
                         cat(sprintf("Error in %s: %s\n", method_name, e$message))
                       })

                       return(method_result)
                     }

  # Stop parallel backend
  stopCluster(cl)

  colnames(results) <- c("prevalence",
                         "anticipated c-stat",
                         "ndev",
                         "method",
                         "calibration slope",
                         "calibration in the large",
                         "auc",
                         "rmspe",
                         paste0("varsel", 1:n.para),
                         "option")

  return(results)
}

