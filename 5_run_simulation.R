source("1_sim_parameters.R")
source("2_data_generation.R")
source("3_performance_measures.R")
source("4_projpred_functions.R")

# Register one parallel backend for cv_varsel(parallel = TRUE)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

proj_results <- foreach(
  i = 1:4,
  .combine = rbind,
  .packages = c("rstanarm", "projpred", "pROC", "dplyr", "mvtnorm", "speedglm"),
  .export = c("run_projpred_one", "run_one_prior", "generate_ss", "measures")
) %dopar% {
  run_projpred_one(i, ndev, n.para, beta0, beta, nval)
}

parallel::stopCluster(cl)
foreach::registerDoSEQ()
