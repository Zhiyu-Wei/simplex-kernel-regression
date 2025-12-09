## ============================================================
##  BLRT for nonlinear component (m1, beta1 = 0)
##  using Nadaraya–Watson Dirichlet kernel smoother (NWD)
##  - Varies the nonlinearity level "a"
## ============================================================

## ---- Load shared NWD functions from GitHub (branch = R) ----
## This file should contain at least:
##   d.kernel_fast()
##   build_K_dirichlet_fast()
##   build_K_dirichlet_fast_UX()
##   CV.Best.h()
source("https://raw.githubusercontent.com/Zhiyu-Wei/simplex-kernel-regression/R/BLRT_NWD.R")

## ---- Packages ----
library(foreach)
library(doParallel)
library(MCMCpack)  # for rdirichlet


## ================== Simulation parameters ==================
a_seq <- seq(0, 0.4, length.out = 10)

M <- 1000    # number of Monte Carlo replications per a
B <- 1000    # number of bootstrap replications
N <- 100
alp <- 0     # Dirichlet kernel parameter for NWD
beta1 <- 0   # true coefficient of z1
type1error <- 0.05

## Matrix to store power / selection summaries
power_result <- matrix(0, nrow = length(a_seq), ncol = 6)
colnames(power_result) <- c("M1vsM2", "M1vsM3", "M2vsM3",
                            "多數選擇", "比例", "sig_beta1")
rownames(power_result) <- paste0("a=", round(a_seq, 2))


## ================== Output folder (auto-created) ==================
## Results will be saved under the current working directory
base_path <- file.path(getwd(), "BLRT_NWD_m1_results")

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

cat("Results will be saved to:", base_path, "\n")


## ================== Parallel setup ==================
nCores <- detectCores() - 8
cl     <- makeCluster(nCores)
registerDoParallel(cl)

cat("使用核心數：", getDoParWorkers(), "\n")

total_start <- Sys.time()


## ================== Outer loop: over different a ==================
for (ai in seq_along(a_seq)) {
  a <- a_seq[ai]
  cat("現在跑 a =", a, "（索引 ai =", ai, "）\n")
  
  ## ---- Inner loop (parallel): M Monte Carlo replications ----
  sim_res <- foreach(
    i = 1:M,
    .combine  = rbind,
    .packages = c("MCMCpack"),
    .export   = c("CV.Best.h", "build_K_dirichlet_fast")
  ) %dopar% {
    
    ## Simple seed in each worker
    set.seed(123 + i)
    
    ## ---------- Data generation ----------
    X  <- rdirichlet(N, c(5, 5, 5))
    x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
    colnames(X) <- c("x1", "x2", "x3")
    z1 <- rnorm(N, 0, 0.01)
    
    f_true <- x1 * x2
    y      <- beta1 * z1 + a * f_true + rnorm(N, 0, 0.01)
    linear.data <- data.frame(X, z1 = z1, y = y)
    n <- nrow(linear.data)
    
    ## ---------- Linear models M1, M2 ----------
    model.1 <- lm(y ~ z1, data = linear.data)
    model.2 <- lm(y ~ z1 + x1 + x2 + x3 - 1, data = linear.data)
    
    lrt_result  <- anova(model.1, model.2, test = "LRT")
    p_M1vsM2    <- lrt_result$`Pr(>Chi)`[2]
    sig_beta1_i <- summary(model.1)$coefficients["z1", 4] < 0.05
    
    ## ---------- Semiparametric model M3 (NWD) ----------
    data2 <- as.matrix(linear.data[, c("x1", "x2", "x3")])
    
    ## Bandwidth selection via CV.Best.h()
    h <- CV.Best.h(data2, y, alp)
    
    ## Nadaraya–Watson smoother S (train | train)
    Sij <- build_K_dirichlet_fast(data2, alp, h)
    Sij <- Sij / rowSums(Sij)
    
    I       <- diag(n)
    Wmatrix <- model.matrix(~ z1 - 1, data = linear.data)
    
    ## Common matrices reused in bootstrap
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
    
    ## ---------- Bootstrap loops ----------
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
    
    ## ---------- Model selection ----------
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
      M1vsM2    = p_M1vsM2,
      M1vsM3    = p_M1vsM3,
      M2vsM3    = p_M2vsM3,
      Select    = selected_model,
      sig_beta1 = sig_beta1_i
    )
  }  # end foreach(i)
  
  ## sim_res is an M x 5 data.frame
  df <- sim_res
  
  ## ---------- Compute power / significance proportions ----------
  power_result[ai, 1] <- mean(df$M1vsM2 < 0.05)
  power_result[ai, 2] <- mean(df$M1vsM3 < 0.05)
  power_result[ai, 3] <- mean(df$M2vsM3 < 0.05)
  
  tab_select <- prop.table(table(df$Select))
  maj_name   <- names(tab_select)[which.max(tab_select)]
  maj_prop   <- max(tab_select)
  
  power_result[ai, 4] <- maj_name
  power_result[ai, 5] <- maj_prop
  power_result[ai, 6] <- mean(df$sig_beta1)
  
  ## ---------- Save raw results for this a ----------
  path_ai <- file.path(
    base_path,
    paste0("power_result_a=", round(a, 2), ".csv")
  )
  write.csv(df, path_ai, row.names = FALSE)
}

stopCluster(cl)

total_end <- Sys.time()
cat("總花費時間：", total_end - total_start, "\n")
print(power_result)

## Also save the summary table
path_power <- file.path(
  base_path,
  "power_result_summary_beta1_0.csv"
)
write.csv(power_result, path_power, row.names = TRUE)


## ================== Read back each a and summarize model selection ==================
path_vec <- file.path(
  base_path,
  paste0("power_result_a=", round(a_seq, 2), ".csv")
)

result_list <- vector("list", length(a_seq))

for (i in seq_along(a_seq)) {
  data_i <- read.csv(path_vec[i])
  tab    <- prop.table(table(data_i$Select))
  
  ## Fill missing categories by 0
  full_tab <- setNames(rep(0, 4), c("M1", "M2", "M3", "None"))
  full_tab[names(tab)] <- tab
  
  result_list[[i]] <- c(
    full_tab,
    beta1 = mean(data_i$sig_beta1)
  )
}

result_df <- do.call(rbind, result_list)
result_df <- as.data.frame(result_df)
result_df$a <- a_seq
result_df  <- result_df[, c("a", "M1", "M2", "M3", "None", "beta1")]

print(result_df, digits = 3, row.names = FALSE)
