require(speedglm)
library(speedglm)
measures <- function(yval, p_val) {
  eta_val <- log(p_val/(1 - p_val))
  # Calibration slope
  fitcal <- speedglm(yval ~ eta_val, family = binomial())
  cal_slope <- as.vector(coef(fitcal)[2])
  # Calibration in the large
  off <- speedglm(yval ~ 1, offset = eta_val, family = binomial())
  cal_large <- as.vector(coef(off))
  # AUC
  cstat <- pROC::roc(response = yval, predictor = as.vector(p_val), levels = c(0, 1), direction = "<")
  auc <- as.vector(cstat$auc)
  # Root mean square prediction error (RMSPE)
  rmspe <- sqrt(mean((p_val - yval)^2))
  return(c(cal_slope, cal_large, auc, rmspe))
}
