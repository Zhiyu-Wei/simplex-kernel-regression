# simplex-kernel-regression

This repository contains the R code accompanying my paper **_Kernel Regression with Simplex Predictors_**.  
It includes implementations of several kernel-based estimators for compositional data, bandwidth selection procedures, simulation experiments, and visualization tools for simplex-valued covariates.

---

## Abstract

In the big data era, some variables with positive values are scaled to sum to one within groups to demonstrate group status or protect individual privacy. The scaled numbers fall in a simplex. Using a linear combination of simplex variables as the conditional mean function in a regression is improper because unexpected constraints are implicitly enforced. The kernel approach is a solution.

In this work, we compare three popular kernel approaches, Aitchison’s logratio-based estimator (1985), the Dirichlet kernel estimator, Bouzebda et al. (2024), and the local linear Dirichlet kernel estimator, Genest and Ouimet (2025), to the semiparametric partial linear framework. A comprehensive simulation study is conducted to compare these methods under diverse data-generating scenarios. The results highlight the strengths and limitations of each approach and, in terms of estimation accuracy, indicate that the local linear estimator consistently outperforms the alternatives, with especially pronounced advantages when the compositional covariates lie close to the simplex boundary. Also, we applied bootstrap likelihood ratio tests to distinguish models with different complexities. Finally, we applied the proposed method to a WHO dataset, treating institutional trust and perceived corruption as regular predictors, and the proportions of education levels across countries as simplex predictors. The analysis showed that the kernel estimator outperforms the linear estimator.

**Keywords:** Dirichlet distribution, nonparametric regression, partial linear regression, simplex

---

# Repository Structure

This repository is organized into three main components to clearly separate reusable functions, figure-generation scripts, and reproducible examples:

R/

Contains all core R functions used in the paper, including:

Kernel estimators (Dirichlet kernel, logratio-based kernel, and local linear Dirichlet kernel)

Bandwidth selection procedures based on LOOCV

Utility functions for simplex visualization and ternary plots
These files provide the fundamental methods used throughout the project.

figures/

Includes:

Scripts used to generate the figures shown in the paper

Example plots that demonstrate how the estimators behave under different settings
These files focus on visualization and graphical diagnostics.

example/

Contains reproducible workflows such as:

Simulation study scripts used to evaluate estimator performance

Code for generating tables included in the manuscript

The real data application using the WHO dataset
These scripts allow readers to reproduce the empirical results in the paper.
