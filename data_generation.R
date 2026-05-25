# generate a dataset (predictors from independent standard normal distributions)
generate_ss <- function(n, n.para, beta0, beta){
  x <- mvtnorm::rmvnorm(n, mean = rep(0, n.para), sigma = diag(n.para))
  eta <- rep(beta0, n) + x%*%beta
  p <- 1/(1+exp(-eta))
  y <- stats::rbinom(n, 1, p)
  data <- data.frame(y,x)}


