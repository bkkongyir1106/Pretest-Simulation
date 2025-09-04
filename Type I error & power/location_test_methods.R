## ---- Load dependencies and user functions --------------------------------
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Type I error & power")

suppressPackageStartupMessages({
  require(foreach)
  require(snow)
  require(doSNOW)
  require(parallel)
})

## ---- Parallel setup helpers ----------------------------------------------
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores()
  cores_use <- max(1, cores - cores_reserve)
  if (Sys.info()[["sysname"]] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doSNOW::registerDoSNOW(cl)  
  } else {
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  cl
}
close_cluster <- function(cl) if (!is.null(cl)) parallel::stopCluster(cl)

## ---- Test functions -------------------------------------------------------
# Bootstrap test for the mean (non-studentized)
one_sample_bootstrap_test <- function(x, mu0 = 0, B = 1000) {
  T_obs <- mean(x) - mu0
  x_centered <- x - mean(x) + mu0
  boot_stats <- replicate(B, mean(sample(x_centered, replace = TRUE)) - mu0)
  mean(abs(boot_stats) >= abs(T_obs))
}

# Sign test for the median
sign_test <- function(x, mu0 = 0) {
  x0 <- x - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(1)  
  signs <- sum(x0 > 0)
  binom.test(signs, n_valid, p = 0.5)$p.value
}

## ---- Parameters -----------------------------------------------------------
Nsim         <- 1e4        
alpha        <- 0.05
B            <- 1e3        
effect_size  <- 0.50       
distributions <- c("Normal", "uniform", "exponential", "LogNormal")
sample_size   <- c(8, 10, 15, 20, 25, 30, 50)

# define distribution parameters
rate <- 1
shape <- 2
meanlog <- 0
sdlog <- 1
## ---- Progress bar + cluster ----------------------------------------------
my_cl <- par_set(cores_reserve = 2)
on.exit(close_cluster(my_cl), add = TRUE)

ntasks <- length(sample_size) * length(distributions)
pb <- txtProgressBar(max = ntasks, style = 3)
on.exit(close(pb), add = TRUE)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## ---- Simulation under H0 and H1 ------------------------------------------
system.time({
  sim_out <- foreach(i = seq_along(sample_size),
                     .packages = c("stats"),
                     .export = c("generate_data",
                                 "one_sample_bootstrap_test",
                                 "sign_test",
                                 "alpha", "B", "Nsim", "effect_size",
                                 "rate","shape", "meanlog", "sdlog",
                                 "distributions", "sample_size")) %:%
    
    foreach(j = seq_along(distributions),
            .options.snow = opts,
            .export = character(0)) %dopar% { #don't re-export in inner loop
              
              n    <- sample_size[i]
              dist <- distributions[j]
              
              ## reproducible RNG per task
              set.seed(12345 + i * 1000 + j, kind = "L'Ecuyer-CMRG")
              
              pval_t_H0  <- pval_boot_H0 <- pval_sign_H0 <- numeric(Nsim)
              pval_t_H1  <- pval_boot_H1 <- pval_sign_H1 <- numeric(Nsim)
              
              for (rep in seq_len(Nsim)) {
                x <- generate_data(n, dist)  
                
                # --- Type I error (H0) ---
                pval_t_H0[rep]      <- t.test(x, mu = 0)$p.value
                pval_boot_H0[rep]   <- one_sample_bootstrap_test(x, mu0 = 0, B = B)
                pval_sign_H0[rep]   <- sign_test(x, mu0 = 0)
                
                
                # --- Power (H1) ---
                x1 <- x + effect_size
                pval_t_H1[rep]      <- t.test(x1, mu = 0)$p.value
                pval_boot_H1[rep]   <- one_sample_bootstrap_test(x1, mu0 = 0, B = B)
                pval_sign_H1[rep]   <- sign_test(x1, mu0 = 0)
              }
              
              list(
                error_t_test         = mean(pval_t_H0 < alpha),
                error_bootstrap_test = mean(pval_boot_H0 < alpha),
                error_sign_test      = mean(pval_sign_H0  < alpha),
                
                power_t_test         = mean(pval_t_H1 < alpha),
                power_bootstrap_test = mean(pval_boot_H1 < alpha),
                power_sign_test      = mean(pval_sign_H1 < alpha)
              )
            }
  
  ## close cluster 
  close_cluster(my_cl)
})

## ---- Collect results into arrays -----------------------------------------
make_arr <- function() array(NA,
                             dim = c(length(sample_size), length(distributions)),
                             dimnames = list(sample_size, distributions))

TypeI_error_t.test         <- make_arr()
TypeI_error_bootstrap.test <- make_arr()
TypeI_error_sign.test      <- make_arr()

power_t.test               <- make_arr()
power_bootstrap.test       <- make_arr()
power_sign.test            <- make_arr()

for (i in seq_along(sample_size)) {
  for (j in seq_along(distributions)) {
    TypeI_error_t.test[i, j]        <- sim_out[[i]][[j]]$error_t_test
    TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
    TypeI_error_sign.test[i, j]      <- sim_out[[i]][[j]]$error_sign_test
    
    power_t.test[i, j]               <- sim_out[[i]][[j]]$power_t_test
    power_bootstrap.test[i, j]       <- sim_out[[i]][[j]]$power_bootstrap_test
    power_sign.test[i, j]            <- sim_out[[i]][[j]]$power_sign_test
  }
}

save(TypeI_error_t.test, TypeI_error_bootstrap.test, TypeI_error_sign.test,
     power_t.test, power_bootstrap.test, power_sign.test, file = "one_sample_location_test.RData")
## ---- print tables ----------------------------------------------
cat("Type I error Rates for t-test\n");         print(TypeI_error_t.test)
cat("Type I error Rates for Bootstrap test\n"); print(TypeI_error_bootstrap.test)
cat("Type I error Rates for Sign test\n");      print(TypeI_error_sign.test)

cat("Power for t-test\n");         print(power_t.test)
cat("Power for Bootstrap test\n"); print(power_bootstrap.test)
cat("Power for Sign test\n");      print(power_sign.test)



## Shared styles
.methods <- c("t-test", "Bootstrap", "Sign")
.cols    <- c("blue",  "black", "red")
.pchs    <- c(16, 17, 18)

## Generic panel plotter (used for both Type I error and Power)
plot_metric_panels <- function(sample_size, distributions, matrices_list,
                               methods = .methods, cols = .cols, pchs = .pchs,
                               main_title = "Metric vs Sample Size",
                               ylab = "Rate", ylim = c(0, 1),
                               alpha_line = NULL, legend_title = "Methods",
                               legend_height = 0.28) {
  
  stopifnot(length(matrices_list) == length(methods),
            length(cols) == length(methods),
            length(pchs) == length(methods))
  
  # Layout: grid of distribution panels + legend row
  nD   <- length(distributions)
  nrow <- ceiling(nD / 2)
  ncol <- min(2, nD)
  
  layout_mat <- rbind(
    matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE),
    rep(nrow * ncol + 1, ncol) # legend row
  )
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  layout(layout_mat, heights = c(rep(1, nrow), legend_height))
  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  # Panels
  for (j in seq_along(distributions)) {
    plot(range(sample_size), ylim, type = "n",
         xlab = "Sample size (n)", ylab = ylab, main = distributions[j])
    if (!is.null(alpha_line))
      abline(h = alpha_line, col = "gray50", lty = 3, lwd = 1.5)
   # grid(col = "lightgray", lty = "dotted")
    
    for (m in seq_along(methods)) {
      y <- matrices_list[[m]][, j]
      lines(sample_size, y, type = "b", lwd = 2, col = cols[m], pch = pchs[m])
    }
  }
  
  # Legend row
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = methods, col = cols, pch = pchs,
         lty = 1, lwd = 2, ncol = 3, bty = "n", title = legend_title, cex = 0.95)
  
  mtext(main_title, outer = TRUE, cex = 1.2, font = 2, line = 0.2)
}

## ---------- Type I error plot ----------
plot_type1_panels <- function(sample_size, distributions, alpha,
                              TypeI_error_t.test, TypeI_error_bootstrap.test, TypeI_error_sign.test) {
  type1_list <- list(TypeI_error_t.test,
                     TypeI_error_bootstrap.test,
                     TypeI_error_sign.test)
  plot_metric_panels(sample_size, distributions, type1_list,
                     methods = .methods, cols = .cols, pchs = .pchs,
                     main_title = "Type I Error vs Sample Size",
                     ylab = "Type I error rate", ylim = c(0, 0.75),
                     alpha_line = alpha, legend_title = "Methods")
}

## ---------- Power plot ----------
plot_power_panels <- function(sample_size, distributions,
                              power_t.test, power_bootstrap.test, power_sign.test,
                              effect_size = NA) {
  power_list <- list(power_t.test,
                     power_bootstrap.test,
                     power_sign.test)
  ttl <- if (is.na(effect_size)) "Power vs Sample Size"
  else sprintf("Power vs Sample Size (effect_size = %.2f)", effect_size)
  
  plot_metric_panels(sample_size, distributions, power_list,
                     methods = .methods, cols = .cols, pchs = .pchs,
                     main_title = ttl, ylab = "Power", ylim = c(0, 1),
                     alpha_line = 0.8, legend_title = "Methods")
}

## =========================
## Draw the two figures
## =========================

# Power
pdf("One_sample_power.pdf", width = 10, height = 8)
plot_power_panels(sample_size, distributions,
                  power_t.test, power_bootstrap.test, power_sign.test,
                  effect_size = effect_size)
dev.off()
# Type I Error
pdf("One_sample_typeI_error.pdf", width = 10, height = 8)
plot_type1_panels(sample_size, distributions, alpha,
                  TypeI_error_t.test, TypeI_error_bootstrap.test, TypeI_error_sign.test)
dev.off()
