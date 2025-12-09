# simplex-kernel-regression

This repository contains the R code accompanying my paper **_Kernel Regression with Simplex Predictors_**.  
It includes implementations of several kernel-based estimators for compositional data, bandwidth selection procedures, simulation experiments, and visualization tools for simplex-valued covariates.

---

## Abstract

In the big data era, some variables with positive values are scaled to sum to one within groups to demonstrate group status or protect individual privacy. The scaled numbers fall in a simplex. Using a linear combination of simplex variables as the conditional mean function in a regression is improper because unexpected constraints are implicitly enforced. The kernel approach is a solution.

In this work, we compare three popular kernel approaches, Aitchison’s logratio-based estimator (1985), the Dirichlet kernel estimator, Bouzebda et al. (2024), and the local linear Dirichlet kernel estimator, Genest and Ouimet (2025), to the semiparametric partial linear framework. A comprehensive simulation study is conducted to compare these methods under diverse data-generating scenarios. The results highlight the strengths and limitations of each approach and, in terms of estimation accuracy, indicate that the local linear estimator consistently outperforms the alternatives, with especially pronounced advantages when the compositional covariates lie close to the simplex boundary. Also, we applied bootstrap likelihood ratio tests to distinguish models with different complexities. Finally, we applied the proposed method to a WHO dataset, treating institutional trust and perceived corruption as regular predictors, and the proportions of education levels across countries as simplex predictors. The analysis showed that the kernel estimator outperforms the linear estimator.

**Keywords:** Dirichlet distribution, nonparametric regression, partial linear regression, simplex

---

## Repository Structure

├── R/ # R functions for kernel estimators & bandwidth selection & visualization tools
├── figures/ # Example visualization scripts (ternary plots, etc.)
├── example/ # Usage examples and reproducible workflows
└── README.md # Project description
