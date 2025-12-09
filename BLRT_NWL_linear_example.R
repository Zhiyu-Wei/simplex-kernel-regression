## ============================================================
##  BLRT for semiparametric model with lognormal kernel (NWL)
##  Linear data-generating process (beta1 = 0)
## ============================================================

## ---- Load shared lognormal-kernel functions from GitHub (branch = R) ----
## This file should contain:
##   - lognormal_kernel_matrix()
##   - lognormal_kernel_reg_cv()
source("https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/R/BLRT_NWL.R")

## ---- Packages used in this simulation ----
library(foreach)
library(doParallel)
library(MCMCpack)  # for rdirichlet


## ================== Simulation parameters ==================
r1_seq <- seq(0, 0.4,   length.out = 10)
r2_seq <- seq(0, 0.274, length.out = 10)
r3_seq <- seq(0, -0.131, length.out = 10)

M <- 1000    # number of Monte Carlo replications per (r1, r2, r3)
B <- 1000    # number of bootstrap replications
N <- 100
alp <- 1     # (kept for symmetry with other codes; not used directly here)
beta1 <- 0
type1error <- 0.05

## Matrix to store power / proportions
power_result <- matrix(0, nrow = length(r1_seq), ncol = 9)
colnames(power_result) <- c("M1vsM2", "M1vsM3", "M2vsM3",
                            "多數選擇", "比例",
                            "β1顯著比例", "r1顯著比例",
                            "r2顯著比例", "r3顯著比例")
rownames(power_result) <- paste0(
  "r1=", round(r1_seq, 2),
  "r2=", round(r2_seq, 2),
  "r3=", round(r3_seq, 2)
)


## ================== Output folder (auto-created) ==================
## All results will be saved under the current working directory
base_path <- file.path(getwd(), "BLRT_NWL_linear_results")

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

cat("Results will be saved to:", base_path, "\n")


## ================== Parallel setup ==================
nCores <- detectCores() - 8
cl     <- makeCluster(nCores)
registerDoParallel(cl)

cat("Number of cores used:", getDoParWorkers(), "\n")

total_start <- Sys.time()


## ================== Outer loop: over (r1, r2, r3) ==================
for (ai in seq_along(r1_seq)) {
  r1 <- r1_seq[ai]
  r2 <- r2_seq[ai]
  r3 <- r3_seq[ai]
  
  cat("Now running combination: ai =", ai,
      " r1 =", r1,
      " r2 =", r2,
      " r3 =", r3, "\n")
  
  ## ---- Inner loop (parallel): M Monte Carlo replications ----
  sim_res <- foreach(
    i = 1:M,
    .combine  = rbind,
    .packages = c("MCMCpack")
  ) %dopar% {
    
    ## Simple seed in each worker
    set.seed(123 + i)
    
    ## ---------- Data generation ----------
    X  <- rdirichlet(N, c(5, 5, 5))
    x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
    z1 <- rnorm(N, 0.5, 1)
    
    y  <- beta1 * z1 + r1 * x1 + r2 * x2 + r3 * x3 + rnorm(N, 0, 0.1)
    linear.data <- data.frame(x1 = x1, x2 = x2, x3 = x3, z1 = z1, y = y)
    n <- nrow(linear.data)
    
    ## ---------- Linear models M1, M2 ----------
    model.1 <- lm(y ~ z1, data = linear.data)
    model.2 <- lm(y ~ z1 + x1 + x2 + x3 - 1, data = linear.data)
    
    anova_m1m2 <- anova(model.1, model.2, test = "LRT")
    p_M1vsM2   <- anova_m1m2$`Pr(>Chi)`[2]
    
    sum_m1    <- summary(model.1)
    sum_m2    <- summary(model.2)
    sig_beta1 <- sum_m1$coefficients["z1", 4] < 0.05
    sig_r1    <- sum_m2$coefficients["x1", 4] < 0.05
    sig_r2    <- sum_m2$coefficients["x2", 4] < 0.05
    sig_r3    <- sum_m2$coefficients["x3", 4] < 0.05
    
    ## ---------- Semiparametric model M3 (lognormal kernel NW) ----------
    data2 <- as.matrix(linear.data[, c("x1", "x2", "x3")])
    
    ## Bandwidth selection via lognormal_kernel_reg_cv()
    h <- lognormal_kernel_reg_cv(data2, y)$h_opt
    
    ## Nadaraya–Watson smoother S = K / rowSums(K)
    Sij <- lognormal_kernel_matrix(data2, h)$K
    Sij <- Sij / rowSums(Sij)
    
    I       <- diag(n)
    Wmatrix <- model.matrix(~ z1 - 1, data = linear.data)
    
    ## Matrices shared between the observed fit and bootstrap fits
    A     <- t(Wmatrix) %*% (I - Sij) %*% Wmatrix
    A_inv <- solve(A)
    
    Betahat <- A_inv %*% t(Wmatrix) %*% (I - Sij) %*% y
    Muhat   <- Sij %*% (y - Wmatrix %*% Betahat)
    pred.y  <- Wmatrix %*% Betahat + Muhat
    RSS_1_M3 <- sum((y - pred.y)^2)
    
    ## ---------- Observed lambda statistics ----------
    lambda_obs_M1M3 <- n * log(deviance(model.1) / RSS_1_M3)
    lambda_obs_M2M3 <- n * log(deviance(model.2) / RSS_1_M3)
    
    res1  <- residuals(model.1)
    pred1 <- fitted(model.1)
    res2  <- residuals(model.2)
    pred2 <- fitted(model.2)
    
    ## ---------- Bootstrap loop ----------
    lambda_boot_M1M3 <- numeric(B)
    lambda_boot_M2M3 <- numeric(B)
    
    for (b in 1:B) {
      ## --- M1 vs M3 ---
      y.star1 <- pred1 + sample(res1, n, replace = TRUE)
      RSS0_1  <- deviance(lm(y.star1 ~ z1, data = linear.data))
      
      Betahat1 <- A_inv %*% t(Wmatrix) %*% (I - Sij) %*% y.star1
      Muhat1   <- Sij %*% (y.star1 - Wmatrix %*% Betahat1)
      pred.y1  <- Wmatrix %*% Betahat1 + Muhat1
      RSS1_1   <- sum((y.star1 - pred.y1)^2)
      
      lambda_boot_M1M3[b] <- n * log(RSS0_1 / RSS1_1)
      
      ## --- M2 vs M3 ---
      y.star2 <- pred2 + sample(res2, n, replace = TRUE)
      RSS0_2  <- deviance(
        lm(y.star2 ~ z1 + x1 + x2 + x3 - 1, data = linear.data)
      )
      
      Betahat2 <- A_inv %*% t(Wmatrix) %*% (I - Sij) %*% y.star2
      Muhat2   <- Sij %*% (y.star2 - Wmatrix %*% Betahat2)
      pred.y2  <- Wmatrix %*% Betahat2 + Muhat2
      RSS1_2   <- sum((y.star2 - pred.y2)^2)
      
      lambda_boot_M2M3[b] <- n * log(RSS0_2 / RSS1_2)
    }
    
    p_M1vsM3 <- (1 + sum(lambda_boot_M1M3 >= lambda_obs_M1M3)) / (B + 1)
    p_M2vsM3 <- (1 + sum(lambda_boot_M2M3 >= lambda_obs_M2M3)) / (B + 1)
    
    ## ---------- Model selection rule ----------
    selected_model <- if (p_M1vsM2 > type1error && p_M1vsM3 > type1error) {
      "M1"
    } else if (p_M1vsM2 < type1error && p_M2vsM3 > type1error) {
      "M2"
    } else if (p_M1vsM3 < type1error && p_M2vsM3 < type1error) {
      "M3"
    } else {
      "None"
    }
    
    ## Return one row per replication
    data.frame(
      M1vsM2 = p_M1vsM2,
      M1vsM3 = p_M1vsM3,
      M2vsM3 = p_M2vsM3,
      Select = selected_model,
      Beta1  = sig_beta1,
      r1     = sig_r1,
      r2     = sig_r2,
      r3     = sig_r3
    )
  }  # end foreach
  
  ## sim_res is an M × 8 data.frame
  df <- sim_res
  
  ## ---------- Compute power / proportions ----------
  power_result[ai, 1] <- mean(df$M1vsM2 < 0.05)
  power_result[ai, 2] <- mean(df$M1vsM3 < 0.05)
  power_result[ai, 3] <- mean(df$M2vsM3 < 0.05)
  
  tab_select <- prop.table(table(df$Select))
  maj_name   <- names(tab_select)[which.max(tab_select)]
  maj_prop   <- max(tab_select)
  
  power_result[ai, 4] <- maj_name
  power_result[ai, 5] <- maj_prop
  power_result[ai, 6] <- mean(df$Beta1)
  power_result[ai, 7] <- mean(df$r1)
  power_result[ai, 8] <- mean(df$r2)
  power_result[ai, 9] <- mean(df$r3)
  
  ## ---------- Save raw results for this (r1, r2, r3) ----------
  path_ai <- file.path(
    base_path,
    paste0("power_result_線性係數=",
           round(r1, 2),
           "beta顯著.csv")
  )
  write.csv(df, path_ai, row.names = FALSE)
}

stopCluster(cl)

## ================== Output power_result ==================
print(power_result)
cat("Total elapsed time:", Sys.time() - total_start, "\n")

path_power <- file.path(
  base_path,
  "power_result線性係數.beta顯著.csv"
)
write.csv(power_result, path_power, fileEncoding = "big5")


## ================== Read back each (r1, r2, r3) and summarize ==================
path_vec <- file.path(
  base_path,
  paste0("power_result_線性係數=",
         round(r1_seq, 2),
         "beta顯著.csv")
)

result_list <- vector("list", length(r1_seq))

for (i in seq_along(r1_seq)) {
  data_i <- read.csv(path_vec[i], fileEncoding = "big5")
  
  ## 4th column is Select (or data_i$Select)
  tab <- prop.table(table(data_i[, 4]))
  
  ## Fill missing categories with 0
  full_tab <- setNames(rep(0, 4), c("M1", "M2", "M3", "None"))
  full_tab[names(tab)] <- tab
  
  ## 5th column is Beta1 (TRUE/FALSE); mean() gives proportion significant
  result_list[[i]] <- c(
    full_tab,
    beta1 = mean(data_i[, 5])
  )
}

result_df <- do.call(rbind, result_list)
result_df <- as.data.frame(result_df)

result_df$r1 <- r1_seq
result_df$r2 <- r2_seq
result_df$r3 <- r3_seq

result_df <- result_df[, c("r1", "r2", "r3",
                           "M1", "M2", "M3", "None",
                           "beta1")]

print(result_df, digits = 3, row.names = FALSE)
