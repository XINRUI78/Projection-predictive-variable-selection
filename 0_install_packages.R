options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c(
  "tidyverse",
  "rstanarm",
  "projpred",
  "loo",
  "glmnet",
  "foreach",
  "doParallel",
  "MASS",
  "pROC",
  "mvtnorm",
  "speedglm",
  "remotes"
)

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if(length(missing) > 0)
  install.packages(missing)
