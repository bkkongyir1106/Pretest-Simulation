# Load dependencies and functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

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

# one sample bootstrap function
one_sample_bootstrap_test <- function(x, n_bootstrap){
  # observe test stats
  observe_test_stat <- sqrt(n) * (mean(x) - 0)/sd(x)
  # bootstrap test stats
  x <- x - mean(x) + mu0
  bootstrap_test_stat <- numeric(n_bootstrap)
  for(b in 1 : n_bootstrap){
    sample1 <- sample(x, n, replace = TRUE)
    bootstrap_test_stat[b] <- sqrt(n) * (mean(sample1) - 0)/sd(sample1)
    bootstrap_p_val <- mean(abs(bootstrap_test_stat) >= abs(observe_test_stat))
  }
  return(bootstrap_p_val)
}


# Parameters
{
  Nsim <- 1e4
  alpha <- 0.05
  n_bootstrap <- 1e3
  effect_size <- 0.5  
  distributions <- c("Normal", "Uniform", "Chi_Square", "Exponential") 
  sample_size <- c(8, 10, 15, 20, 25, 30, 50)
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
# One-sample case simulation
# ----------------------------

# Simulation under H0 and H1
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow = opts) %:%
    foreach(dist = distributions) %dopar% {
      set.seed(12345)
      pval_t_H0 <- pval_wilcox_H0 <- pval_boot_H0 <- numeric(Nsim)
      pval_t_H1 <- pval_wilcox_H1 <- pval_boot_H1 <- numeric(Nsim)
      
      for (i in 1:Nsim) {
        x <- generate_data(n, dist)
        
        # Type I error
        pval_t_H0[i] <- t.test(x)$p.value
        pval_wilcox_H0[i] <- wilcox.test(x)$p.value
        pval_boot_H0[i] <- one_sample_bootstrap_test(x, n_bootstrap = n_bootstrap)
        
        # power
        pval_t_H1[i] <- t.test(x + effect_size)$p.value
        pval_wilcox_H1[i] <- wilcox.test(x + effect_size)$p.value
        pval_boot_H1[i] <- one_sample_bootstrap_test(x + effect_size, n_bootstrap = n_bootstrap)
      }
      
      list(
        error_t_test = mean(pval_t_H0 < alpha),
        error_wilcox_test = mean(pval_wilcox_H0 < alpha),
        error_bootstrap_test = mean(pval_boot_H0 < alpha),
        power_t_test = mean(pval_t_H1 < alpha),
        power_wilcox_test = mean(pval_wilcox_H1 < alpha),
        power_bootstrap_test = mean(pval_boot_H1 < alpha)
      )
    }
  
 # close_cluster(my_cl)
  
  # Store results
  errorvec <- numeric(length(sample_size) * length(distributions))
  TypeI_error_t.test <- TypeI_error_wilcox.test <- TypeI_error_bootstrap.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  power_t.test <- power_wilcox.test <- power_bootstrap.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      TypeI_error_t.test[i, j] <- sim_out[[i]][[j]]$error_t_test
      TypeI_error_wilcox.test[i, j] <- sim_out[[i]][[j]]$error_wilcox_test
      TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
      power_t.test[i, j] <- sim_out[[i]][[j]]$power_t_test
      power_wilcox.test[i, j] <- sim_out[[i]][[j]]$power_wilcox_test
      power_bootstrap.test[i, j] <- sim_out[[i]][[j]]$power_bootstrap_test
    }
  }
})

# Print results
cat("Type I error Rates for t-test\n"); print(TypeI_error_t.test)
cat("Type I error Rates for Wilcoxon-test\n"); print(TypeI_error_wilcox.test)
cat("Type I error Rates for Bootstrap test\n"); print(TypeI_error_bootstrap.test)
cat("Power for t-test\n"); print(power_t.test)
cat("Power for Wilcoxon-test\n"); print(power_wilcox.test)
cat("Power for Bootstrap test\n"); print(power_bootstrap.test)

# Save results
save(Nsim, n_bootstrap, effect_size, distributions, sample_size, 
     TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test, 
     power_t.test, power_wilcox.test, power_bootstrap.test, 
     file = "one_sample_test.RData")


# Create plots for Type I error rates
# -------------------------------------------
create_plots_typeI <- function() {
  # Layout: 2 rows x 2 columns for plots, 1 row for shared legend
  layout_matrix <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE)
  layout(mat = layout_matrix, heights = c(1, 1, 0.3))
  
  # Plotting setup
  par(mar = c(4, 4, 2, 1))
  colors <- c("blue", "green", "red")
  pchs <- c(16, 17, 21)
  tests <- c("t-test", "Wilcoxon test", "Bootstrap test")
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get data for each test
    y_t <- TypeI_error_t.test[, j]
    y_w <- TypeI_error_wilcox.test[, j]
    y_b <- TypeI_error_bootstrap.test[, j]
    
    # Determine ylim
    y_max <- max(y_t, y_w, y_b, 0.2, na.rm = TRUE)
    
    # Create base plot
    plot(sample_size, y_t, type = "b", col = colors[1], lty = 1, pch = pchs[1], lwd = 2, ylim = c(0, y_max),xlab = "Sample Size", ylab = "Type I Error Rate", main = dist_name)
    
    # Add other tests
    lines(sample_size, y_w, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_b, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
    
    # Add reference line at alpha = 0.05
    abline(h = alpha, col = "gray", lty = 3, lwd = 2)
  }
  
  # Legend (bottom row)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = tests, title = "Test Method", col = colors, lty = 1, lwd = 2, pch = pchs, horiz = FALSE, ncol = 3, bty = "b", cex = 1.0)
}

# Create plots for Power
create_plots_power <- function() {
  # Layout: 2 rows x 2 columns for plots, 1 row for shared legend
  layout_matrix <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE)
  layout(mat = layout_matrix, heights = c(1, 1, 0.3))
  
  # Plotting setup
  par(mar = c(4, 4, 2, 1))
  colors <- c("blue", "green", "red")
  pchs <- c(16, 17, 21, 15)
  tests <- c("t-test", "Wilcoxon test", "Bootstrap test")
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get data for each test
    y_t <- power_t.test[, j]
    y_w <- power_wilcox.test[, j]
    y_b <- power_bootstrap.test[, j]
    
    # Determine ylim
    y_max <- 1
    
    # Create base plot
    plot(sample_size, y_t, type = "b", col = colors[1], lty = 1, pch = pchs[1], lwd = 2, ylim = c(0, y_max), xlab = "Sample Size", ylab = "Power", main = dist_name)
    
    # Add other tests
    lines(sample_size, y_w, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_b, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
  }
  
  # Legend (bottom row)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = tests, title = "Test Method", col = colors, lty = 1, lwd = 2, pch = pchs, horiz = FALSE, ncol = 3, bty = "b", cex = 1.0)
}

# Create PDF files with the plots
pdf("one_sample_TypeI_error_rates.pdf", width = 10, height = 8)
create_plots_typeI()
dev.off()

pdf("one_sample_Power_rates.pdf", width = 10, height = 8)
create_plots_power()
dev.off()

## -------------------------------------------------------
##                Two-sample case simulation
## ------------------------------------------------------

# Bootstrap two-sample test
bootstrap_two_sample_test <- function(x, y, effect_size = 0,  n_bootstrap = 1000) {
  observe_test_stat <- (mean(x) - mean(y) -effect_size)/sqrt(var(x)/length(x) + var(y)/length(y))
  combined <- c(x, y)
  n1 <- length(x)
  n2 <- length(y)
  
  # Center under H0: mean(x) - mean(y) = effect_size
  # Set group means to (grand_mean + effect_size/2) and (grand_mean - effect_size/2)
  grand_mean <- mean(combined)
  x_centered <- x - mean(x) + (grand_mean + effect_size/2)
  y_centered <- y - mean(y) + (grand_mean - effect_size/2)
  
  boot_stats <- replicate(n_bootstrap, {
    x_star <- sample(x_centered, n1, replace = TRUE)
    y_star <- sample(y_centered, n2, replace = TRUE)
    (mean(x_star) - mean(y_star) - effect_size)/sqrt(var(x_star)/length(x_star) + var(y_star)/length(y_star))
  })
  
  p_val <- mean(abs(boot_stats) >= abs(observe_test_stat))
  return(p_val)
}

# Simulation under H0 and H1
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow = opts) %:%
    foreach(dist = distributions) %dopar% {
      set.seed(12345)
      pval_t_H0 <- pval_wilcox_H0 <- pval_boot_H0 <- numeric(Nsim)
      pval_t_H1 <- pval_wilcox_H1 <- pval_boot_H1 <- numeric(Nsim)
      
      for (i in 1:Nsim) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        
        # Type I error
        pval_t_H0[i] <- t.test(x, y)$p.value
        pval_wilcox_H0[i] <- wilcox.test(x, y)$p.value
        pval_boot_H0[i] <- bootstrap_two_sample_test(x, y, effect_size = 0, n_bootstrap = n_bootstrap)
        
        # power
        pval_t_H1[i] <- t.test(x, y + effect_size)$p.value
        pval_wilcox_H1[i] <- wilcox.test(x, y + effect_size)$p.value
        pval_boot_H1[i] <- bootstrap_two_sample_test(x, y, effect_size = effect_size, n_bootstrap = n_bootstrap)
      }
      
      list(
        error_t_test = mean(pval_t_H0 < alpha),
        error_wilcox_test = mean(pval_wilcox_H0 < alpha),
        error_bootstrap_test = mean(pval_boot_H0 < alpha),
        power_t_test = mean(pval_t_H1 < alpha),
        power_wilcox_test = mean(pval_wilcox_H1 < alpha),
        power_bootstrap_test = mean(pval_boot_H1 < alpha)
      )
    }
  
  close_cluster(my_cl)
  
  # Store results
  TypeI_error_t.test <- TypeI_error_wilcox.test <- TypeI_error_bootstrap.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  power_t.test <- power_wilcox.test <- power_bootstrap.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      TypeI_error_t.test[i, j] <- sim_out[[i]][[j]]$error_t_test
      TypeI_error_wilcox.test[i, j] <- sim_out[[i]][[j]]$error_wilcox_test
      TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
      power_t.test[i, j] <- sim_out[[i]][[j]]$power_t_test
      power_wilcox.test[i, j] <- sim_out[[i]][[j]]$power_wilcox_test
      power_bootstrap.test[i, j] <- sim_out[[i]][[j]]$power_bootstrap_test
    }
  }
})

# Print results
cat("Type I error Rates for t-test\n"); print(TypeI_error_t.test)
cat("Type I error Rates for Wilcoxon-test\n"); print(TypeI_error_wilcox.test)
cat("Type I error Rates for Bootstrap test\n"); print(TypeI_error_bootstrap.test)
cat("Power for t-test\n"); print(power_t.test)
cat("Power for Wilcoxon-test\n"); print(power_wilcox.test)
cat("Power for Bootstrap test\n"); print(power_bootstrap.test)

# Save results
save(Nsim, B, effect_size, distributions, sample_size, 
     TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test, 
     power_t.test, power_wilcox.test, power_bootstrap.test, 
     file = "Two_sample_test.RData")

# Create plots for two-sample case
pdf("two_sample_TypeI_error_rates.pdf", width = 10, height = 8)
create_plots_typeI()
dev.off()

pdf("two_sample_Power_rates.pdf", width = 10, height = 8)
create_plots_power()
dev.off()

