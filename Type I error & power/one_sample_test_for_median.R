## ---- Setup & deps -------------------------------------------------------
# source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

library(MASS)
library(foreach)
library(doSNOW)   # uniform backend w/ progress
library(doRNG)    # reproducible parallel RNG
library(parallel)

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Type I error & power")

## ---- Parallel helpers ---------------------------------------------------
par_set <- function(cores_reserve = 2) {
  total <- max(1L, parallel::detectCores(logical = TRUE))
  use   <- max(1L, total - as.integer(cores_reserve))
  # PSOCK works cross-platform and supports doSNOW progress callbacks
  cl <- parallel::makePSOCKcluster(use)
  doSNOW::registerDoSNOW(cl)
  # Load packages + source scripts on workers
  parallel::clusterEvalQ(cl, {
    library(MASS); library(VGAM); library(LaplacesDemon)
  })
  # If your sourced files define functions needed on workers, source them there too:
  # parallel::clusterEvalQ(cl, {
  #   source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")
  #   source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
  # })
  cl
}
close_cluster <- function(cl) {
  try(parallel::stopCluster(cl), silent = TRUE)
  foreach::registerDoSEQ()
}

## ---- Tests (unchanged logic) -------------------------------------------
sign_test <- function(x, mu0 = 0) {
  x0 <- x - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  list(p.value = stats::binom.test(signs, n_valid, p = 0.5)$p.value)
}

one_sample_bootstrap_t <- function(x, n_bootstrap, mu0 = 0){
  n <- length(x)
  obs <- stats::median(x) - mu0
  x_null <- x - stats::median(x) + mu0
  tstar <- numeric(n_bootstrap)
  for(b in seq_len(n_bootstrap)){
    xb <- sample(x_null, n, replace = TRUE)
    mb <- stats::median(xb)
    M <- 200L
    se_b <- stats::sd(replicate(M, stats::median(sample(xb, n, replace = TRUE))))
    tstar[b] <- (mb - mu0) / se_b
  }
  M0 <- 1000L
  se_obs <- stats::sd(replicate(M0, stats::median(sample(x, n, replace = TRUE))))
  tobs <- obs / se_obs
  mean(abs(tstar) >= abs(tobs))
}

one_sample_perm <- function(x, mu0 = 0, B = 5000) {
  n <- length(x)
  obs <- stats::median(x)
  perm <- replicate(B, {
    s <- sample(c(-1, 1), n, replace = TRUE)
    x_star <- mu0 + s * (x - mu0)
    stats::median(x_star)
  })
  mean(abs(perm) >= abs(obs))
}

boxcox_simple <- function(x) {
  shift <- 0
  if (min(x) <= 0) {
    shift <- abs(min(x)) + 0.001
    x <- x + shift
  }
  bc <- suppressMessages(MASS::boxcox(x ~ 1, lambda = seq(-2, 2, 0.1), plotit = FALSE))
  lambda <- bc$x[which.max(bc$y)]
  transformed_x <- if (lambda == 0) log(x) else (x^lambda - 1) / lambda
  list(transformed_x = transformed_x, lambda = lambda, shift = shift)
}

bc_t_test <- function(x, mu0 = 0) {
  bc <- boxcox_simple(x)
  median_shifted <- mu0 + bc$shift
  transformed_median <- if (bc$lambda == 0) log(median_shifted) else (median_shifted^bc$lambda - 1) / bc$lambda
  pval <- stats::t.test(bc$transformed_x, mu = transformed_median)$p.value
  list(p.value = pval)
}

## ---- Parameters ---------------------------------------------------------
Nsim         <- 1e2
alpha        <- 0.05
n_bootstrap  <- 1e2
effect_size  <- 0.5
distributions <- c("normal", "laplace")
sample_size   <- c(10, 20, 30, 40, 50)

## ---- Progress + cluster -------------------------------------------------
my_cl <- par_set(cores_reserve = 2)
# unified cleanup (cluster, progress bar, backend)
pb <- txtProgressBar(max = length(sample_size) * length(distributions), style = 3)
progress_count <- local({ i <- 0L; function(...) { i <<- i + 1L; setTxtProgressBar(pb, i) } })
opts <- list(progress = progress_count)

on.exit({
  try(close(pb), silent = TRUE)
  try(doParallel::stopImplicitCluster(), silent = TRUE)
  close_cluster(my_cl)
}, add = FALSE)

## ---- Simulation ---------------------------------------------------------
# Make RNG reproducible *and* independent across workers
set.seed(12345)
registerDoRNG(12345)

system.time({
  # Export locally defined helpers to workers
  .exports_vec <- c("sign_test","one_sample_bootstrap_t","one_sample_perm",
                    "boxcox_simple","bc_t_test","alpha","Nsim","n_bootstrap",
                    "effect_size","distributions","sample_size")
  # If generate_data() lives in sourced files that are *not* sourced on workers,
  # add "generate_data" to .exports_vec (or source those files via clusterEvalQ above).
  .exports_vec <- unique(c(.exports_vec, "generate_data"))
  
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon","VGAM","MASS"),
                     .options.snow = opts,
                     .export = .exports_vec,
                     .inorder = TRUE) %dopar% {
                       foreach(dist = distributions,
                               .packages = c("LaplacesDemon","VGAM","MASS"),
                               .export = .exports_vec,
                               .inorder = TRUE) %do% {
                                 
                                 pval_t_H0 <- pval_sign_H0 <- pval_boot_H0 <- pval_perm_H0 <- pval_bc_H0 <- numeric(Nsim)
                                 pval_t_H1 <- pval_sign_H1 <- pval_boot_H1 <- pval_perm_H1 <- pval_bc_H1 <- numeric(Nsim)
                                 
                                 for (i in seq_len(Nsim)) {
                                   x <- generate_data(n, dist)
                                   
                                   # H0
                                   pval_sign_H0[i] <- sign_test(x, mu0 = 0)$p.value
                                   pval_t_H0[i]    <- stats::t.test(x, mu = 0)$p.value
                                   pval_boot_H0[i] <- one_sample_bootstrap_t(x, n_bootstrap = n_bootstrap, mu0 = 0)
                                   pval_perm_H0[i] <- one_sample_perm(x, mu0 = 0, B = n_bootstrap)
                                   pval_bc_H0[i]   <- bc_t_test(x, mu0 = 0)$p.value
                                   
                                   # H1
                                   x_shifted <- x + effect_size
                                   pval_sign_H1[i] <- sign_test(x_shifted, mu0 = 0)$p.value
                                   pval_t_H1[i]    <- stats::t.test(x_shifted, mu = 0)$p.value
                                   pval_boot_H1[i] <- one_sample_bootstrap_t(x_shifted, n_bootstrap = n_bootstrap, mu0 = 0)
                                   pval_perm_H1[i] <- one_sample_perm(x_shifted, mu0 = 0, B = n_bootstrap)
                                   pval_bc_H1[i]   <- bc_t_test(x_shifted, mu0 = 0)$p.value
                                 }
                                 
                                 list(
                                   error_t_test         = mean(pval_t_H0    < alpha),
                                   error_sign_test      = mean(pval_sign_H0 < alpha),
                                   error_bootstrap_test = mean(pval_boot_H0 < alpha),
                                   error_perm_test      = mean(pval_perm_H0 < alpha),
                                   error_bc_test        = mean(pval_bc_H0   < alpha),
                                   power_sign_test      = mean(pval_sign_H1 < alpha),
                                   power_t_test         = mean(pval_t_H1    < alpha),
                                   power_bootstrap_test = mean(pval_boot_H1 < alpha),
                                   power_perm_test      = mean(pval_perm_H1 < alpha),
                                   power_bc_test        = mean(pval_bc_H1   < alpha)
                                 )
                               }
                     }
  
  # Collect results
  errorvec <- numeric(length(sample_size) * length(distributions))
  TypeI_error_t.test        <- array(errorvec, dim = c(length(sample_size), length(distributions)),
                                     dimnames = list(sample_size, distributions))
  TypeI_error_sign.test     <- TypeI_error_t.test
  TypeI_error_bootstrap.test<- TypeI_error_t.test
  TypeI_error_perm.test     <- TypeI_error_t.test
  TypeI_error_bc.test       <- TypeI_error_t.test
  
  power_t.test              <- TypeI_error_t.test
  power_sign.test           <- TypeI_error_t.test
  power_bootstrap.test      <- TypeI_error_t.test
  power_perm.test           <- TypeI_error_t.test
  power_bc.test             <- TypeI_error_t.test
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      TypeI_error_t.test[i, j]         <- sim_out[[i]][[j]]$error_t_test
      TypeI_error_sign.test[i, j]      <- sim_out[[i]][[j]]$error_sign_test
      TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
      TypeI_error_perm.test[i, j]      <- sim_out[[i]][[j]]$error_perm_test
      TypeI_error_bc.test[i, j]        <- sim_out[[i]][[j]]$error_bc_test
      
      power_sign.test[i, j]            <- sim_out[[i]][[j]]$power_sign_test
      power_t.test[i, j]               <- sim_out[[i]][[j]]$power_t_test
      power_bootstrap.test[i, j]       <- sim_out[[i]][[j]]$power_bootstrap_test
      power_perm.test[i, j]            <- sim_out[[i]][[j]]$power_perm_test
      power_bc.test[i, j]              <- sim_out[[i]][[j]]$power_bc_test
    }
  }
})

## ---- Print & Save -------------------------------------------------------
cat("Type I error Rates for Sign test\n"); print(TypeI_error_sign.test)
cat("Type I error Rates for t test\n"); print(TypeI_error_t.test)
cat("Type I error Rates for Bootstrap test\n"); print(TypeI_error_bootstrap.test)
cat("Type I error Rates for Permutation test\n"); print(TypeI_error_perm.test)
cat("Type I error Rates for Box-Cox t-test\n"); print(TypeI_error_bc.test)

cat("Power for Sign test\n"); print(power_sign.test)
cat("Power for t test\n"); print(power_t.test)
cat("Power for Bootstrap test\n"); print(power_bootstrap.test)
cat("Power for Permutation test\n"); print(power_perm.test)
cat("Power for Box-Cox t-test\n"); print(power_bc.test)

save(Nsim, n_bootstrap, effect_size, distributions, sample_size,
     TypeI_error_t.test, TypeI_error_sign.test, TypeI_error_bootstrap.test,
     TypeI_error_perm.test, TypeI_error_bc.test,
     power_t.test, power_sign.test, power_bootstrap.test,
     power_perm.test, power_bc.test,
     file = "one_sample_median_tests.RData")

## ---- Plotting (unchanged) ----------------------------------------------
create_simple_combined_plot <- function() {
  par(mfrow = c(2, 2), oma = c(0.5, 0, 2, 0))
  colors <- c("blue", "green", "red", "purple", "brown")
  pchs <- c(16, 17, 15, 18, 21)
  tests <- c("t-test","Sign test", "Bootstrap test", "Permutation test", "Box-Cox t-test")
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    y_t   <- TypeI_error_t.test[, j]
    y_sign<- TypeI_error_sign.test[, j]
    y_boot<- TypeI_error_bootstrap.test[, j]
    y_perm<- TypeI_error_perm.test[, j]
    y_bc  <- TypeI_error_bc.test[, j]
    y_max <- max(y_t, y_sign, y_boot, y_perm, y_bc, 0.08, na.rm = TRUE)
    
    plot(sample_size, y_t, type = "b", col = colors[1], lty = 1, pch = pchs[1], lwd = 2,
         ylim = c(0, y_max), xlab = "Sample Size", ylab = "Type I Error Rate",
         main = paste(dist_name, "- Type I Error"))
    lines(sample_size, y_sign, type = "b", col = colors[2], pch = pchs[2], lwd = 2)
    lines(sample_size, y_boot, type = "b", col = colors[3], pch = pchs[3], lwd = 2)
    lines(sample_size, y_perm, type = "b", col = colors[4], pch = pchs[4], lwd = 2)
    lines(sample_size, y_bc,   type = "b", col = colors[5], pch = pchs[5], lwd = 2)
    abline(h = alpha, col = "gray", lty = 3, lwd = 2)
  }
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    y_t   <- power_t.test[, j]
    y_sign<- power_sign.test[, j]
    y_boot<- power_bootstrap.test[, j]
    y_perm<- power_perm.test[, j]
    y_bc  <- power_bc.test[, j]
    
    plot(sample_size, y_t, type = "b", col = colors[1], lty = 1, pch = pchs[1], lwd = 2,
         ylim = c(0, 1), xlab = "Sample Size", ylab = "Power",
         main = paste(dist_name, "- Power"))
    lines(sample_size, y_sign, type = "b", col = colors[2], pch = pchs[2], lwd = 2)
    lines(sample_size, y_boot, type = "b", col = colors[3], pch = pchs[3], lwd = 2)
    lines(sample_size, y_perm, type = "b", col = colors[4], pch = pchs[4], lwd = 2)
    lines(sample_size, y_bc,   type = "b", col = colors[5], pch = pchs[5], lwd = 2)
  }
  
  title("Comparison of Type I Error Rates and Power for One Sample Median Tests (test 2)",
        outer = TRUE, cex.main = 1.2)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = tests, col = colors, lty = 1, lwd = 2,
         pch = pchs, horiz = TRUE, bty = "n", cex = 0.9, xpd = TRUE)
}

pdf("one_sample_median_simple_combined.pdf", width = 9, height = 8)
create_simple_combined_plot()
dev.off()
