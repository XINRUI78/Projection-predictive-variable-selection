# Register one parallel backend for cv_varsel(parallel = TRUE)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

proj_results <- do.call(
  rbind,
  lapply(1:10, function(i) {
    run_projpred_one(i, ndev, n.para, beta0, beta, nval)
  })
)

parallel::stopCluster(cl)
foreach::registerDoSEQ()
