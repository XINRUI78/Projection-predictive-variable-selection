source("0_install_packages.R")
source("1_sim_parameters.R")
source("2_data_generation.R")
source("3_performance_measures.R")
source("4_nor.R")
#source("4_projpred_functions.R")

cl <- parallel::makeCluster(30)
doParallel::registerDoParallel(cl)

proj_results <- foreach(
  i = 1:100,
  .combine = rbind,
  .packages = c("rstanarm", "projpred", "pROC", "mvtnorm", "speedglm"),
  .export = c("run_one_prior", "generate_ss", "measures")
) %dopar% {
  #run_projpred_one(i, ndev, n.para, beta0, beta, nval)
  run_one_prior(i, ndev, nval, n.para, beta0, beta)
  }

parallel::stopCluster(cl)
foreach::registerDoSEQ()

colnames(proj_results) <- c(
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

# Save output
write.csv(proj_results, "proj_results.csv", row.names = FALSE)
saveRDS(proj_results, "proj_results.rds")
