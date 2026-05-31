# Set simulation parameters
n.para <- 30
prev <- 0.3
c <- 0.8
nval <- 10000
percentage <- c(0.1, 0.2, 0.2, 0.5)

# Define number of predictors with relative strengths
strong <- percentage[1] * n.para  # 10 strong predictive variables
medium <- percentage[2] * n.para  # 30% medium predictive variables
weak <- percentage[3] * n.para    # 30% weak predictive variables
noise <- percentage[4] * n.para   # 30% noise predictive variables

# Assign relative strengths
weights <- c(rep(1, strong), rep(0.5, medium), rep(0.25, weak), rep(0, noise))

# Calculate recommended sample size
#library(samplesizedev)

#rss <- samplesizedev(outcome = "Binary", S = 0.9, phi = prev, c = c, p = n.para)
#ndev <- rss$sim # recommended sample size
ndev <- 1538
ndev1 <- round(ndev/2) # half the recommended sample size
ndev2 <- round(ndev/4) # one-quarter recommended sample size

# The optimizer function that obtains the true coefficients
opt_beta <- function(n.para, prev, c, weights) {
  # Generate predictors (X) from multivariate normal distribution
  n = 500000
  x <- rmvnorm(n, mean = rep(0, n.para), sigma = diag(n.para))
  objective <- function(para){
    beta0 <- para[1]  # Intercept
    s <- para[2]      # Scaling factor
    beta1 <- s * weights
    eta <- rep(beta0, n) + x %*% beta1
    p <- 1/(1+exp(-eta))
    y <- rbinom(n, 1, p)
    pest <- mean(y)
    cstat <- roc(response = as.vector(y), predictor = as.vector(p), levels = c(0, 1), direction = "<")
    cest <- as.vector(cstat$auc)
    return((pest - prev)^2 + (cest - c)^2)
  }
  # Initial guesses for beta0 and s
  initial_para <- c(-2, 1)
  tol = 1e-6
  # Perform optimization
  result <- optim(
    par = initial_para,
    fn = objective,
    method = "Nelder-Mead",
    control = list(abstol = tol)
  )

  # Extract optimized coefficients
  beta0_opt <- result$par[1]
  s_opt <- result$par[2]
  beta1_opt <- s_opt * weights

  list(
    beta0 = beta0_opt,
    beta1 = beta1_opt,
    s = s_opt
  )
}

# Obtain the true coefficients
opt_beta <- opt_beta(n.para, prev, c, weights)
beta0 <- opt_beta$beta0
beta <- opt_beta$beta1
