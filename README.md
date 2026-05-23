# Projection Predictive Variable Selection Simulation Study

This repository contains R code for a simulation study evaluating Bayesian projection predictive variable selection for logistic regression models.

The study compares Bayesian reference models with different priors and evaluates projected submodels using external validation datasets.

---

# Overview

For each simulation replicate:

1. Development dataset and validation dataset are generated from a predefined logistic regression data-generating mechanism.
2. A full Bayesian logistic regression reference model is fitted.
3. Projection predictive variable selection is performed using forward search with 5-fold cross validation.
4. Projected submodels are evaluated on the validation dataset.

The procedure is repeated across 100 simulated datasets.

---

# Statistical Framework

## Reference Models

Bayesian logistic regression reference models are fitted using:

```r
stan_glm()
```

from the `rstanarm` package.

Two prior distributions are considered:

### Normal prior

\[
\beta_j \sim \mathcal{N}(0,1)
\]

### Laplace prior

\[
\beta_j \sim \text{Laplace}(0,1)
\]

Hamiltonian Monte Carlo sampling is used to obtain posterior draws from the reference model.

---

# Projection Predictive Variable Selection

Projection predictive variable selection is implemented using:

```r
cv_varsel()
```

from the `projpred` package.

The procedure:

1. Fits a full Bayesian reference model.
2. Draws posterior samples from the reference model.
3. Projects posterior draws onto candidate submodels.
4. Minimizes the Kullback–Leibler divergence between:
   - the predictive distribution of the reference model
   - the predictive distribution of the submodel.
5. Evaluates predictive performance using cross-validation.

Forward search is used to construct submodels sequentially.

---

# Mathematical Description

For each posterior draw of the reference model:

\[
\theta_*^{(s)}
\]

the reference model produces predictive probabilities:

\[
p(\tilde y \mid \theta_*^{(s)})
\]

For each candidate submodel, projected parameters:

\[
\theta_\perp^{(s)}
\]

are obtained by minimizing:

\[
\text{KL}\left(
p(\tilde y \mid \theta_*^{(s)})
\;\|\;
p(\tilde y \mid \theta_\perp^{(s)})
\right)
\]

The projected parameter draws:

\[
\{
\theta_\perp^{(1)},
\dots,
\theta_\perp^{(S)}
\}
\]

form the projected posterior distribution.

Predictions are averaged across projected posterior draws using:

```r
proj_linpred(..., integrated = TRUE)
```

---

# Cross-Validation

The simulation uses:

```r
validate_search = TRUE
```

which repeats the entire variable-selection search inside each cross-validation fold.

Current settings:

- 5-fold cross-validation
- forward search
- parallel computation

This produces more reliable estimates of predictive performance but substantially increases computation time.

---

# Model Selection

Two projected submodels are selected for each reference model:

## Suggested model

Obtained using:

```r
suggest_size()
```

This selects the smallest model whose predictive performance is sufficiently close to the reference model.

## Best model

Defined as the model with the maximum expected log predictive density (ELPD).

---

# Performance Measures

Models are evaluated on the external validation dataset using:

- Calibration slope
- Calibration-in-the-large
- Area under the ROC curve (AUC)
- Root mean squared prediction error (RMSPE)

Variable-selection performance is also recorded using binary indicators showing whether each predictor was selected.

---

# Main Function

## `run_projpred_one()`

Runs one simulation replicate.

### Inputs

| Argument | Description |
|---|---|
| `i` | Simulation seed/index |
| `ndev` | Development sample size |
| `n.para` | Number of candidate predictors |
| `beta0` | True intercept |
| `beta` | True coefficient vector |
| `nval` | Validation sample size |

### Output

Returns a matrix containing results for:

| Method |
|---|
| Normal prior reference model |
| Normal projected suggested-size model |
| Normal projected best-size model |
| Laplace prior reference model |
| Laplace projected suggested-size model |
| Laplace projected best-size model |

---

# Running the Simulation

Example:

```r
proj_results <- do.call(
  rbind,
  lapply(1:100, function(i) {
    run_projpred_one(i, ndev, n.para, beta0, beta, nval)
  })
)
```

---

# Parallel Computing

Parallel computation is implemented using:

```r
foreach
doParallel
```

Parallelization is used during:

- fold-wise variable selection
- performance evaluation

inside `cv_varsel()`.

---

# Required Packages

```r
library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)
```

Additional helper functions may require:

```r
library(pROC)
library(mvtnorm)
```

---

# Saving Results

## Save as CSV

```r
write.csv(
  proj_results,
  "projpred_results.csv",
  row.names = FALSE
)
```

## Save as R object

```r
saveRDS(
  proj_results,
  "projpred_results.rds"
)
```

---

# Notes

- `validate_search = TRUE` is computationally intensive.
- Parallel computation may substantially increase memory usage.
- Running Stan chains and projpred parallelization simultaneously may exhaust system memory.
- If memory issues occur:
  - reduce the number of workers,
  - reduce the number of chains,
  - reduce `nterms_max`,
  - or set `parallel = FALSE`.

---

# Author

Simulation study for Bayesian projection predictive variable selection in logistic regression.
