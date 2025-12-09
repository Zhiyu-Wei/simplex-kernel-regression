## ============================================================
##  BLRT for nonlinear component (m1, beta1 = 0)
##  using local linear Dirichlet kernel smoother
##  - Varies the nonlinearity level "a"
## ============================================================

## ---- Load shared functions from GitHub (branch = R) ----
## This file should contain:
##   d.kernel_fast()
##   build_K_dirichlet_fast()
##   build_K_dirichlet_fast_UX()
##   WLS_best_h()
##   build_S_locallinear_dirichlet()
source("https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/R/BLRT_LLD.R")

## ---- Packages used in the main simulation ----
library(foreach)
library(doParallel)
library(MCMCpack)  # for rdirichlet


## ============================================================
##  Simulation settings
## ============================================================
a_seq      <- seq(0, 0.4, length.out = 10)   # strength of nonlinear term
M          <- 1000    # number of Monte Carlo replications per a
B          <- 1000    # number of bootstrap replications for BLRT
N          <- 100     # sample size
alp        <- 1       # Dirichlet kernel parameter alpha
beta1      <- 0       # true coefficient of z1
type1error <- 0.05    # nominal significance level


## ============================================================
##  Parallel setup
## ============================================================
nCores <- detectCores() - 8    # leave some cores free
cl     <- makeCluster(nCores)
registerDoParallel(cl)

cat("Number of cores used:", getDoParWorkers(), "\n")


## ============================================================
##  Storage for power / type I error summary
## ============================================================
power_result <- matrix(0, nrow = length(a_seq), ncol = 6)
colnames(power_result) <- c("M1vsM2", "M1vsM3", "M2vsM3",
                            "多數選擇", "比例", "sig_beta1")
rownames(power_result) <- paste0("a=", round(a_seq, 2))


## ============================================================
##  Output folder (created automatically in working directory)
## ============================================================
base_path <- file.path(getwd(), "BLRT_LLD_nonlinear_results")

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

cat("Results will be saved to:", base_path, "\n")

total_start <- Sys.time()


## ============================================================
##  Outer loop over different values of a
## ============================================================
for (ai in seq_along(a_seq)) {
  a <- a_seq[ai]
  cat("Now running a =", a, "(index ai =", ai, ")\n")
  
  ## ----------------------------------------------------------
  ##  Inner loop (parallel): M Monte Carlo replications
  ## ----------------------------------------------------------
  sim_res <- foreach(
    i = 1:M,
    .combine  = rbind,
    .packages = c("MCMCpack")   # add more packages here if needed
  ) %dopar% {
    
    ## Simple random seed (not perfect for strict reproducibility under parallelism,
    ## but usually sufficient for simulation work)
    set.seed(123 + i)
    
    ## ------------------------------------------------------
    ## 1. Data generation
    ## ------------------------------------------------------
    X  <- rdirichlet(N, c(5, 5, 5))
    x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
    colnames(X) <- c("x1", "x2", "x3")
    
    ## z1 is almost constant (small variance), so its effect
    ## should be hard to detect under H0: beta1 = 0
    z1 <- rnorm(N, 0, 0.01)
    
    ## Nonlinear function in x: f(x) = x1 * x2
    f_true <- x1 * x2
    
    ## Data-generating model: y = beta1*z1 + a*f(x) + error
    y <- beta1 * z1 + a * f_true + rnorm(N, 0, 0.01)
    
    linear.data <- data.frame(X, z1 = z1, y = y)
    n <- nrow(linear.data)
    
    ## ------------------------------------------------------
    ## 2. Linear models M1 and M2
    ## ------------------------------------------------------
    ## M1: y ~ z1
    model.1 <- lm(y ~ z1, data = linear.data)
    ## M2: y ~ z1 + x1 + x2 + x3 (no intercept)
    model.2 <- lm(y ~ z1 + x1 + x2 + x3 - 1, data = linear.data)
    
    ## Likelihood ratio test: M1 vs M2
    lrt_result <- anova(model.1, model.2, test = "LRT")
    p_M1M2     <- lrt_result$`Pr(>Chi)`[2]
    
    ## Significance of beta1 in M1
    sig_b1 <- summary(model.1)$coefficients["z1", 4] < 0.05
    
    ## ------------------------------------------------------
    ## 3. Semiparametric model M3 (local linear Dirichlet)
    ##    Using x = (x1, x2, x3) as kernel variables
    ##    and (x1, x2) as local-linear covariates.
    ## ------------------------------------------------------
    data2 <- as.matrix(linear.data[, c("x1", "x2", "x3")])
    
    ## Bandwidth selection using WLS_best_h()
    h <- WLS_best_h(data2, y, alpha = 1)
    
    X_kern_train <- data2                        # kernel part (x1, x2, x3)
    X_loc_train  <- data2[, 1:2, drop = FALSE]   # local-linear part (x1, x2)
    
    ## Smoother matrix S(y) for nonparametric component
    Sij <- build_S_locallinear_dirichlet(
      X_kern = X_kern_train,
      X_loc  = X_loc_train,
      alpha  = alp,
      h      = h
    )
    
    I       <- diag(n)
    Wmatrix <- model.matrix(~ z1 - 1, data = linear.data)
    
    ## Common matrices (used in both observed fit and bootstrap loop)
    A     <- t(Wmatrix) %*% (I - Sij) %*% Wmatrix
    A_inv <- solve(A)
    
    Betahat <- A_inv %*% t(Wmatrix) %*% (I - Sij) %*% y
    Muhat   <- Sij %*% (y - Wmatrix %*% Betahat)
    pred.y  <- Wmatrix %*% Betahat + Muhat
    RSS_1_M3 <- sum((y - pred.y)^2)
    
    ## ------------------------------------------------------
    ## 4. Observed test statistics lambda(M1 vs M3), lambda(M2 vs M3)
    ## ------------------------------------------------------
    lambda_obs_M1M3 <- n * log(deviance(model.1) / RSS_1_M3)
    lambda_obs_M2M3 <- n * log(deviance(model.2) / RSS_1_M3)
    
    res.1  <- residuals(model.1)
    pred.1 <- fitted(model.1)
    res.2  <- residuals(model.2)
    pred.2 <- fitted(model.2)
    
    ## ------------------------------------------------------
    ## 5. Bootstrap loop (within each worker)
    ## ------------------------------------------------------
    lambda_boot_M1M3 <- numeric(B)
    lambda_boot_M2M3 <- numeric(B)
    
    for (b in 1:B) {
      ## --- Bootstrap for M1 vs M3 ---
      y.star1 <- pred.1 + sample(res.1, n, replace = TRUE)
      RSS_0_1 <- deviance(lm(y.star1 ~ z1, data = linear.data))
      
      Betahat1 <- A_inv %*% t(Wmatrix) %*% (I - Sij) %*% y.star1
      Muhat1   <- Sij %*% (y.star1 - Wmatrix %*% Betahat1)
      pred.y1  <- Wmatrix %*% Betahat1 + Muhat1
      RSS_1_1  <- sum((y.star1 - pred.y1)^2)
      
      lambda_boot_M1M3[b] <- n * log(RSS_0_1 / RSS_1_1)
      
      ## --- Bootstrap for M2 vs M3 ---
      y.star2 <- pred.2 + sample(res.2, n, replace = TRUE)
      RSS_0_2 <- deviance(
        lm(y.star2 ~ z1 + x1 + x2 + x3 - 1, data = linear.data)
      )
      
      Betahat2 <- A_inv %*% t(Wmatrix) %*% (I - Sij) %*% y.star2
      Muhat2   <- Sij %*% (y.star2 - Wmatrix %*% Betahat2)
      pred.y2  <- Wmatrix %*% Betahat2 + Muhat2
      RSS_1_2  <- sum((y.star2 - pred.y2)^2)
      
      lambda_boot_M2M3[b] <- n * log(RSS_0_2 / RSS_1_2)
    }
    
    ## ------------------------------------------------------
    ## 6. Bootstrap p-values
    ## ------------------------------------------------------
    p_M1M3 <- (1 + sum(lambda_boot_M1M3 >= lambda_obs_M1M3)) / (B + 1)
    p_M2M3 <- (1 + sum(lambda_boot_M2M3 >= lambda_obs_M2M3)) / (B + 1)
    
    ## ------------------------------------------------------
    ## 7. Model selection rule
    ## ------------------------------------------------------
    select <- if (p_M1M2 > type1error && p_M1M3 > type1error) {
      "M1"
    } else if (p_M1M2 < type1error && p_M2M3 > type1error) {
      "M2"
    } else if (p_M1M3 < type1error && p_M2M3 < type1error) {
      "M3"
    } else {
      "None"
    }
    
    ## Return one row per replication
    data.frame(
      M1vsM2    = p_M1M2,
      M1vsM3    = p_M1M3,
      M2vsM3    = p_M2M3,
      Select    = select,
      sig_beta1 = as.numeric(sig_b1)
    )
  }  # end foreach
  
  
  ## ========================================================
  ##  Aggregate results for this value of a
  ## ========================================================
  result.M1vsM2 <- sim_res$M1vsM2
  result.M1vsM3 <- sim_res$M1vsM3
  result.M2vsM3 <- sim_res$M2vsM3
  model.select  <- sim_res$Select
  sig_beta1     <- sim_res$sig_beta1
  
  ## Save raw results for this a as CSV
  final <- data.frame(
    M1vsM2    = result.M1vsM2,
    M1vsM3    = result.M1vsM3,
    M2vsM3    = result.M2vsM3,
    Select    = model.select,
    sig_beta1 = sig_beta1
  )
  
  path_a <- file.path(
    base_path,
    paste0("power_result_a=", round(a, 2), ".csv")
  )
  write.csv(final, path_a, row.names = FALSE)
  
  ## Compute power / selection proportion / beta1 significance
  power_result[ai, 1] <- mean(result.M1vsM2 < 0.05)
  power_result[ai, 2] <- mean(result.M1vsM3 < 0.05)
  power_result[ai, 3] <- mean(result.M2vsM3 < 0.05)
  
  tab_select <- prop.table(table(model.select))
  max_idx    <- which.max(tab_select)
  power_result[ai, 4] <- names(tab_select)[max_idx]
  power_result[ai, 5] <- as.numeric(tab_select[max_idx])
  power_result[ai, 6] <- mean(sig_beta1)
}

stopCluster(cl)

total_end <- Sys.time()
cat("Total elapsed time:", total_end - total_start, "\n")
print(power_result)


## ============================================================
##  Read back each a result and summarize selection proportions
## ============================================================
path_vec <- file.path(
  base_path,
  paste0("power_result_a=", round(a_seq, 2), ".csv")
)

result_list <- vector("list", length(a_seq))

for (i in seq_along(a_seq)) {
  data_i <- read.csv(path_vec[i])
  
  ## Column "Select" stores the selected model for each replication
  tab <- prop.table(table(data_i$Select))
  
  ## Ensure all four categories exist in the vector (fill missing with zero)
  full_tab <- setNames(rep(0, 4), c("M1", "M2", "M3", "None"))
  full_tab[names(tab)] <- tab
  
  ## Proportion of times beta1 was significant
  result_list[[i]] <- c(full_tab,
                        beta1 = mean(data_i$sig_beta1))
}

result_df <- do.call(rbind, result_list)
result_df <- as.data.frame(result_df)
result_df$a <- a_seq
result_df  <- result_df[, c("a", "M1", "M2", "M3", "None", "beta1")]

print(result_df, digits = 3, row.names = FALSE)

