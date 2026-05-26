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

The main performance measures include:

- Calibration slope
- Calibration-in-the-large
- AUC 
- RMSPE
- Number of selected predictors
- Number of true predictors selected

---



For each simulated dataset, the code:

1. Generates a development dataset.
2. Generates an independent validation dataset.
3. Fits prediction models using different variable-selection methods.
4. Evaluates model performance on the validation dataset.
5. Saves performance measures and selected predictors.



---

## Main R Files

### `1_sim_parameters.R`

This file defines the simulation settings and generates the true regression coefficients used in the simulation study.

The purpose of this script is to:

1. Define the simulation scenario (30 predictors, prevalence = 0.3, expected AUC = 0.8, strength of strong, medium, weak, noise predictors = (1, 0.5, 0.25, 0), the percentage of corresponding type = (0.1, 0.2, 0.2, 0.5))
2. Calculate the required development sample size (ndev)
3. Generate the true logistic regression coefficients that achieve the desired:
   - outcome prevalence
   - AUC

### `2_data_generation.R`

This function is used to generate development datasets and validation datasets (size of 10000) for the simulation study. The predictors are from independent standard normal distributions.

### `3_performance_measures.R`

This function calculates validation performance measures:
- calibration slope
- calibration-in-the-large
- AUC
- RMSPE
using the predicted probabilities and observed outcomes from validation datasets.

### `4_projpred_functions.R`

Contains functions for Bayesian projection predictive variable selection.

Main functions include:

```r
run_one_prior()
run_projpred_one()
```
#function run_one_prior()
The code fits Bayesian logistic regression reference models using all variables and a specified prior distribution. It then applies projection predictive variable selection using forward search with 5-fold cross-validation.
It is designed so the same code can be reused for different priors.

#function run_projpred_one()
It runs one complete simulation replicate.
For one simulated development and validation dataset, it fits and evaluates:
1. Normal-prior reference model
2. Normal-prior projected suggested model (model size is selected using 1-se rule)
3. Normal-prior projected best model (model size is selected using the best elpd)
4. Laplace-prior reference model
5. Laplace-prior projected suggested model (model size is selected using 1-se rule)
6. Laplace-prior projected best model (model size is selected using the best elpd)


### `5_run_simulation.R`

This is the main script used to run the simulation 100 times.

### `frequentist_methods.R`

Contains code for traditional frequentist methods, such as:

- full maximum likelihood logistic regression
- univariable screening
- backward elimination
- LASSO
- modified LASSO
- ridge regression
- random forest, if included

Run this file after loading the data-generation and performance-measure functions.

