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
