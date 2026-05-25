# Register one parallel backend for cv_varsel(parallel = TRUE)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

proj_results <- foreach(
  i = 1:100,
  .combine = rbind,
  .packages = c("rstanarm", "projpred", "pROC", "dplyr"),
  .export = c("run_projpred_one", "run_one_prior", "generate_ss", "measures")
) %dopar% {
  run_projpred_one(i, ndev, n.para, beta0, beta, nval)
}

parallel::stopCluster(cl)
foreach::registerDoSEQ()
