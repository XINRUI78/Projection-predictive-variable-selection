# Projection Predictive Variable Selection Simulation Study

This repository contains R code for a simulation study comparing Bayesian projection predictive variable selection using different reference model priors.

## Overview

The simulation evaluates Bayesian logistic regression reference models followed by projection predictive variable selection. For each simulated dataset, the code fits:

1. A full Bayesian reference model with a Normal prior.
2. Projected submodels selected from the Normal-prior reference model.
3. A full Bayesian reference model with a Laplace prior.
4. Projected submodels selected from the Laplace-prior reference model.

The workflow is repeated over multiple simulated development and validation datasets.

## Main Function

The main function is:

```r
run_projpred_one(i, ndev, n.para, beta0, beta, nval)
