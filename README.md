# simplex-kernel-regression

This repository is organized into three main components to separate reusable functions, visualization scripts, and reproducible examples.  
Several scripts and development files are maintained in additional branches (such as **R** and **figures**).  
If a specific file is not present in the main branch, please check these branches for the latest version.

---

## Abstract

In the big data era, some variables with positive values are scaled to sum to one within groups to demonstrate group status or protect individual privacy. The scaled numbers fall in a simplex. Using a linear combination of simplex variables as the conditional mean function in a regression is improper because unexpected constraints are implicitly enforced. The kernel approach is a solution.

In this work, we compare three popular kernel approaches, Aitchison’s logratio-based estimator (1985), the Dirichlet kernel estimator, Bouzebda et al. (2024), and the local linear Dirichlet kernel estimator, Genest and Ouimet (2025), to the semiparametric partial linear framework. A comprehensive simulation study is conducted to compare these methods under diverse data-generating scenarios. The results highlight the strengths and limitations of each approach and, in terms of estimation accuracy, indicate that the local linear estimator consistently outperforms the alternatives, with especially pronounced advantages when the compositional covariates lie close to the simplex boundary. Also, we applied bootstrap likelihood ratio tests to distinguish models with different complexities. Finally, we applied the proposed method to a WHO dataset, treating institutional trust and perceived corruption as regular predictors, and the proportions of education levels across countries as simplex predictors. The analysis showed that the kernel estimator outperforms the linear estimator.

**Keywords:** Dirichlet distribution, nonparametric regression, partial linear regression, simplex

---

## Repository Structure

This repository is organized into three main components to separate reusable functions, visualization scripts, and reproducible examples.

### **R** [👉 Click here: R code](https://github.com/Zhiyu-Wei/simplex-kernel-regression/tree/R)

This folder contains all core R functions used in the project, including:

- Kernel estimation methods (Dirichlet kernel, logratio-based estimator, and local linear Dirichlet kernel)
- Bandwidth selection procedures based on leave-one-out cross-validation (LOOCV)
- Visualization and simplex-related utilities for ternary plotting

These files provide the essential methodological tools used across the analyses.

### **figures** [👉 Click here: figures](https://github.com/Zhiyu-Wei/simplex-kernel-regression/tree/figure)
This folder includes:

- Scripts used to generate figures referenced in the paper
- Example visualization outputs such as ternary plots and smoothing illustrations

Files in this directory focus on graphical diagnostics and visual demonstrations.

### **example** [👉 Click here: Simulation Results & Application](https://github.com/Zhiyu-Wei/simplex-kernel-regression/tree/example)
This folder contains reproducible workflows and applied scripts, such as:

- Simulation study code evaluating estimator performance
- Scripts for producing tables included in the manuscript
- The real data application using the WHO dataset

These scripts allow readers to reproduce all empirical results from the study.

