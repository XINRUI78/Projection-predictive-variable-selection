source("0_install_packages.R")
source("1_sim_parameters.R")
source("2_data_generation.R")
source("3_performance_measures.R")
source("4_nor.R")

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
run_dir <- Sys.getenv("RUN_DIR")

proj_results <- run_nor(
  i = task_id,
  ndev = ndev,
  nval = nval,
  n.para = n.para,
  beta0 = beta0,
  beta = beta
)

write.csv(
  proj_results,
  sprintf("array_nor_%03d.csv", task_id),
  row.names = FALSE
)
