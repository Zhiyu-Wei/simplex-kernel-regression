# Simulation Results & Application

## Overview

This simulation study evaluates how different smoothing methods behave in a semi-parametric regression setting and how well the Bootstrap Likelihood Ratio Test (BLRT) can select the correct model.  
We focus on three main questions:

---

### (1) Effect of the centering parameter \(a\)
[👉 Click here: Effect of a & Comparison](https://github.com/Zhiyu-Wei/simplex-kernel-regression/tree/example/Effect-of-a-%26-Comparison)


We first investigate how the centering parameter \(a\) affects semi-parametric estimation when using Dirichlet-based kernel smoothers.

All simulations follow a unified setup that includes four data-generating mechanisms (two linear and two nonlinear), each examined under two types of simplex-valued covariate designs:

- **Center-concentrated design:** most observations lie around the middle of the simplex  
- **Boundary-concentrated design:** many observations are located near the edges

Because \(a\) controls the centering and localization of the Dirichlet kernel, different values of \(a\) may influence both the nonlinear smoothing behavior and the estimation of the linear component in the partial linear model.

For each scenario, we evaluate the effect of \(a\) for both NWD and LLD, focusing on:

- the accuracy of the nonlinear estimator, measured by the integrated squared error (ISE), and  
- the accuracy and stability of the estimated linear coefficient.

This analysis highlights how the choice of \(a\) affects Dirichlet-type smoothers and whether its influence differs between NWD and LLD.

---

### (2) Comparison of smoothing methods

Using the same simulation settings, we compare three smoothing approaches:

- **NWD** (Dirichlet Nadaraya–Watson)  
- **NWL** (Lognormal Nadaraya–Watson)  
- **LLD** (Local Linear Dirichlet)

For each scenario and replication, we fit the semi-parametric model with each method and evaluate:

- the estimation of the linear component, and  
- the accuracy of the estimated nonlinear function, measured by MISE.

This comparison shows how the different smoothing methods behave under both linear and nonlinear structures, and how stable they remain when the covariates concentrate near the simplex boundary.

---

### (3) Model selection performance of the BLRT

Finally, we evaluate how well the BLRT distinguishes between competing models across different data-generating conditions.

We vary two types of signal strength:

1. the linear part, $ \mathbf z_i^{\top}\boldsymbol\beta $, and  
2. the nonlinear part, $ \mu(\mathbf x_i) $.

By adjusting these signal levels, we observe how the BLRT responds as model complexity changes and how effectively it identifies the correct specification.

---
