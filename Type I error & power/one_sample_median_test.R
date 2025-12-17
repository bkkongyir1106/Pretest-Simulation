# Load dependencies and functions
#source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
library(MASS)

# kill_workers <- function() {
#   try(doParallel::stopImplicitCluster(), TRUE)
#   try({
#     cl <- parallel::getDefaultCluster()
#     if (!is.null(cl)) parallel::stopCluster(cl)
#   }, TRUE)
#   foreach::registerDoSEQ()
# }

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Type I error & power")

# Parallel setup
par_set <- function(cores_reserve = 2) {
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  return(cl)
}
close_cluster <- function(cl) parallel::stopCluster(cl)

# test functions

## ---------- Sign test for the median ---------
sign_test <- function(x, mu0 = 0) {
  x0 <- x - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  return(list(p.value = binom.test(signs, n_valid, p = 0.5)$p.value))
}

# # one sample bootstrap function
# one_sample_bootstrap_test <- function(x, n_bootstrap, mu0 = 0){
#   n <- length(x)
#   # observe test stats
#   observe_test_stat <- median(x)
#   # bootstrap test stats
#   x <- x - median(x) + mu0
#   bootstrap_test_stat <- numeric(n_bootstrap)
#   for(b in 1 : n_bootstrap){
#     sample1 <- sample(x, n, replace = TRUE)
#     bootstrap_test_stat[b] <- median(sample1)
#   }
#   bootstrap_p_val <- mean(abs(bootstrap_test_stat) >= abs(observe_test_stat))
#   return(bootstrap_p_val)
# }


one_sample_bootstrap_t <- function(x, n_bootstrap, mu0 = 0){
  n <- length(x)
  obs <- median(x) - mu0
  
  # null-resampled data
  x_null <- x - median(x) + mu0
  
  # outer bootstrap under H0: get T* = (med* - mu0)/SE*(med*)
  tstar <- numeric(n_bootstrap)
  for(b in 1:n_bootstrap){
    xb <- sample(x_null, n, replace = TRUE)
    mb <- median(xb)
    
    # inner SE estimate via quick bootstrap (small M works fine)
    M <- 200L
    se_b <- sd(replicate(M, median(sample(xb, n, replace = TRUE))))
    tstar[b] <- (mb - mu0) / se_b
  }
  
  # observed t using SE from the data itself
  M0 <- 1000L
  se_obs <- sd(replicate(M0, median(sample(x, n, replace = TRUE))))
  tobs <- obs / se_obs
  
  mean(abs(tstar) >= abs(tobs))
}


# -----------One Sample Permutation test ------
one_sample_perm <- function(x, mu0 = 0, B = 5000) {
  n <- length(x)
  obs <- median(x)
  perm <- replicate(B, {
    s <- sample(c(-1, 1), n, replace = TRUE)
    x_star <- mu0 + s * (x - mu0)
    median(x_star)
  })
  mean(abs(perm) >= abs(obs))
}

# ----------------------------------------------------
# Box-Cox transformation 
boxcox_simple <- function(x) {
  # Shift x to be positive
  shift <- 0
  if (min(x) <= 0) {
    shift <- abs(min(x)) + 0.001
    x <- x + shift
  }
  
  # Find optimal lambda 
  bc <- suppressMessages(boxcox(x ~ 1, lambda = seq(-2, 2, 0.1), plotit = FALSE))
  lambda <- bc$x[which.max(bc$y)]
  
  # Transform x
  if (lambda == 0) {
    transformed_x <- log(x)
  } else {
    transformed_x <- (x^lambda - 1) / lambda
  }
  
  return(list(
    transformed_x = transformed_x,
    lambda = lambda,
    shift = shift
  ))
}

# ---------------------------------------------------
# T-test on median using Box-Cox transformed x
bc_t_test <- function(x, mu0 = 0) {
  # Transform x
  bc <- boxcox_simple(x)
  
  # Transform the median under H0
  median_shifted <- mu0 + bc$shift
  if (bc$lambda == 0) {
    transformed_median <- log(median_shifted)
  } else {
    transformed_median <- (median_shifted^bc$lambda - 1) / bc$lambda
  }
  
  # T-test on transformed x against transformed median
  pval <- t.test(bc$transformed_x, mu = transformed_median)$p.value
  return(list(p.value = pval))
}

# Parameters
{
  Nsim <- 1e3
  alpha <- 0.05
  n_bootstrap <- 1e3
  effect_size <- 0.5  
  distributions <- c("normal", "laplace") 
  sample_size <- c(10, 20,  30, 40, 50)
}

# Progress bar
{
  my_cl <- par_set(cores_reserve = 2)
  on.exit(close_cluster(my_cl))  # ensures clean exit
  ntasks <- length(sample_size) * length(distributions)
  pb <- txtProgressBar(max = ntasks, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
}

# ----------------------------
# One-sample case simulation for median testing
# ----------------------------

on.exit({
  try(doParallel::stopImplicitCluster(), TRUE)
  try(parallel::stopCluster(cl), TRUE)
  foreach::registerDoSEQ()
}, add = TRUE)


# Simulation under H0 and H1
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM", "MASS"),
                     .options.snow = opts) %:%
    foreach(dist = distributions) %dopar% {
      set.seed(12345)
      pval_t_H0 <- pval_sign_H0 <- pval_boot_H0 <- pval_perm_H0 <- pval_bc_H0 <- numeric(Nsim)
      pval_t_H1 <- pval_sign_H1 <- pval_boot_H1 <- pval_perm_H1 <- pval_bc_H1 <- numeric(Nsim)
      
      for (i in 1:Nsim) {
        # Generate data centered around median 0
        x <- generate_data(n, dist)
        
        # Type I error (testing H0: median = 0)
        pval_sign_H0[i] <- sign_test(x, mu0 = 0)$p.value
        pval_t_H0[i] <- t.test(x, mu0 = 0)$p.value
        pval_boot_H0[i] <- one_sample_bootstrap_t(x, n_bootstrap = n_bootstrap, mu0 = 0)
        pval_perm_H0[i] <- one_sample_perm(x, mu0 = 0, B = n_bootstrap)
        pval_bc_H0[i] <- bc_t_test(x, mu0 = 0)$p.value
        
        # Power (testing H0: median = 0 when true median = effect_size)
        x_shifted <- x + effect_size
        pval_sign_H1[i] <- sign_test(x_shifted, mu0 = 0)$p.value
        pval_t_H1[i] <- t.test(x_shifted, mu0 = 0)$p.value
        pval_boot_H1[i] <- one_sample_bootstrap_t(x_shifted, n_bootstrap = n_bootstrap, mu0 = 0)
        pval_perm_H1[i] <- one_sample_perm(x_shifted, mu0 = 0, B = n_bootstrap)
        pval_bc_H1[i] <- bc_t_test(x_shifted, mu0 = 0)$p.value
      }
      
      list(
        error_t_test = mean(pval_t_H0 < alpha),
        error_sign_test = mean(pval_sign_H0 < alpha),
        error_bootstrap_test = mean(pval_boot_H0 < alpha),
        error_perm_test = mean(pval_perm_H0 < alpha),
        error_bc_test = mean(pval_bc_H0 < alpha),
        power_sign_test = mean(pval_sign_H1 < alpha),
        power_t_test = mean(pval_t_H1 < alpha),
        power_bootstrap_test = mean(pval_boot_H1 < alpha),
        power_perm_test = mean(pval_perm_H1 < alpha),
        power_bc_test = mean(pval_bc_H1 < alpha)
      )
    }
  
  # Store results
  errorvec <- numeric(length(sample_size) * length(distributions))
  TypeI_error_t.test <- TypeI_error_sign.test <- TypeI_error_bootstrap.test <- TypeI_error_perm.test <- TypeI_error_bc.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  power_t.test <-power_sign.test <- power_bootstrap.test <- power_perm.test <- power_bc.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      TypeI_error_t.test[i, j] <- sim_out[[i]][[j]]$error_t_test
      TypeI_error_sign.test[i, j] <- sim_out[[i]][[j]]$error_sign_test
      TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
      TypeI_error_perm.test[i, j] <- sim_out[[i]][[j]]$error_perm_test
      TypeI_error_bc.test[i, j] <- sim_out[[i]][[j]]$error_bc_test
      power_sign.test[i, j] <- sim_out[[i]][[j]]$power_sign_test
      power_t.test[i, j] <- sim_out[[i]][[j]]$power_t_test
      power_bootstrap.test[i, j] <- sim_out[[i]][[j]]$power_bootstrap_test
      power_perm.test[i, j] <- sim_out[[i]][[j]]$power_perm_test
      power_bc.test[i, j] <- sim_out[[i]][[j]]$power_bc_test
    }
  }
})

# close cluster
close_cluster(my_cl)

# Print results
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

# Save results
save(Nsim, n_bootstrap, effect_size, distributions, sample_size, TypeI_error_t.test,
     TypeI_error_sign.test, TypeI_error_bootstrap.test, TypeI_error_perm.test, TypeI_error_bc.test,
     power_t.test, power_sign.test, power_bootstrap.test, power_perm.test, power_bc.test, 
     file = "one_sample_median_tests.RData")

# ---------------------------------------------------------
# Create plots for Type I error rates
create_simple_combined_plot <- function() {
  # Set up 2x2 layout with space for title
  par(mfrow = c(2, 2), oma = c(0.5, 0, 2, 0))  # Add outer margin at top for title
  
  # Plotting setup
  colors <- c("blue", "green", "red", "purple", "brown")
  pchs <- c(16, 17, 15, 18, 21)
  tests <- c("t-test","Sign test", "Bootstrap test", "Permutation test", "Box-Cox t-test")
  
  # Type I Error plots (top row)
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get Type I error data for each test
    y_t <- TypeI_error_t.test[, j]
    y_sign <- TypeI_error_sign.test[, j]
    y_boot <- TypeI_error_bootstrap.test[, j]
    y_perm <- TypeI_error_perm.test[, j]
    y_bc <- TypeI_error_bc.test[, j]
    
    # Determine ylim
    y_max <- max(y_t, y_sign, y_boot, y_perm, y_bc, 0.08, na.rm = TRUE)
    
    # Create base plot
    plot(sample_size, y_t, type = "b", col = colors[1], lty = 1, pch = pchs[1], lwd = 2, 
         ylim = c(0, y_max), xlab = "Sample Size", ylab = "Type I Error Rate", 
         main = paste(dist_name, "- Type I Error"))
    
    # Add other tests
    lines(sample_size, y_sign, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_boot, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
    lines(sample_size, y_perm, type = "b", col = colors[4], lty = 1, pch = pchs[4], lwd = 2)
    lines(sample_size, y_bc, type = "b", col = colors[5], lty = 1, pch = pchs[5], lwd = 2)
    
    # Add reference line at alpha = 0.05
    abline(h = alpha, col = "gray", lty = 3, lwd = 2)
  }
  
  # Power plots (bottom row)
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get power data for each test
    y_t <- power_t.test[, j]
    y_sign <- power_sign.test[, j]
    y_boot <- power_bootstrap.test[, j]
    y_perm <- power_perm.test[, j]
    y_bc <- power_bc.test[, j]
    
    # Create base plot
    plot(sample_size, y_t, type = "b", col = colors[1], lty = 1, pch = pchs[1], lwd = 2, 
         ylim = c(0, 1), xlab = "Sample Size", ylab = "Power", 
         main = paste(dist_name, "- Power"))
    
    # Add other tests
    lines(sample_size, y_sign, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_boot, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
    lines(sample_size, y_perm, type = "b", col = colors[4], lty = 1, pch = pchs[4], lwd = 2)
    lines(sample_size, y_bc, type = "b", col = colors[5], lty = 1, pch = pchs[5], lwd = 2)
  }
  
  # Add main title
  title("Comparison of Type I Error Rates and Power for One Sample Median Tests(test 2)", outer = TRUE, cex.main = 1.2)
  
  # Add single legend to the entire plot
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = tests, col = colors, lty = 1, lwd = 2,  
         pch = pchs, horiz = TRUE, bty = "n", cex = 0.9, xpd = TRUE)
}

# Create PDF file with the simple combined plot
pdf("one_sample_median_simple_combined.pdf", width = 9, height = 8)
create_simple_combined_plot()
dev.off()
