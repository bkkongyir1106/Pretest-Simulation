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

# Bootstrap test for the mean
one_sample_bootstrap_test <- function(x, mu0 = 0, B = 1000) {
  T_obs <- mean(x) - mu0
  x_centered <- x - mean(x) + mu0
  boot_stats <- replicate(B, mean(sample(x_centered, replace = TRUE)) - mu0)
  p_val <- mean(abs(boot_stats) >= abs(T_obs))
  return(p_val)
}

# Parameters
{
  Nsim <- 1e4
  alpha <- 0.05
  B <- 1e3
  effect_size <- 0.5  
  distributions <- c("Normal", "Uniform", "t", "Exponential", "Chi-Square", "LogNormal")
  sample_size <- c(8, 10, 15, 20, 25, 30, 50)
}

# Define plot function
create_plots <- function(results, y_var, title, ylab, filename, 
                         colors = c("blue", "green", "red"), 
                         pchs = c(16, 17, 21),
                         tests = c("t-test", "Wilcoxon test", "Bootstrap test"),
                         hline = NULL) {
  
  pdf(filename, width = 10, height = 8)
  
  # Layout: 2 rows x 3 columns for plots, 1 row for shared legend
  layout_matrix <- matrix(1:9, nrow = 3, byrow = TRUE)
  layout(mat = layout_matrix, heights = c(1, 1, 0.3))
  
  # Plotting setup
  par(mar = c(4, 4, 2, 1))
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get data for each test
    y_t <- results[[1]][, j]
    y_w <- results[[2]][, j]
    y_b <- results[[3]][, j]
    
    # Determine ylim
    y_max <- if(y_var == "power") 1 else max(y_t, y_w, y_b, 0.2)
    
    # Create base plot
    plot(sample_size, y_t, 
         type = "b", 
         col = colors[1], 
         lty = 1, 
         pch = pchs[1],
         lwd = 2,
         ylim = c(0, y_max),
         xlab = "Sample Size", 
         ylab = ylab,
         main = dist_name)
    
    # Add other tests
    lines(sample_size, y_w, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_b, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
    
    # Add reference line if specified
    if(!is.null(hline)) abline(h = hline, col = "gray", lty = 3, lwd = 2)
  }
  
  # Legend (bottom row)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = tests,
         title = "Test Method", 
         col = colors, 
         lty = 1,
         lwd = 2,
         pch = pchs,
         horiz = TRUE,
         bty = "n", 
         cex = 1.2)
  
  dev.off()
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
        pval_boot_H0[i] <- one_sample_bootstrap_test(x, mu0 = 0, B)
        
        # power
        pval_t_H1[i] <- t.test(x + effect_size)$p.value
        pval_wilcox_H1[i] <- wilcox.test(x + effect_size)$p.value
        pval_boot_H1[i] <- one_sample_bootstrap_test(x + effect_size, mu0 = 0, B)
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
save(Nsim, B, effect_size, distributions, sample_size, 
     TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test, 
     power_t.test, power_wilcox.test, power_bootstrap.test, 
     file = "one_sample_test.RData")

# Create plots for one-sample case
create_plots(list(TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test),
             y_var = "error", title = "Type I Error Rates", ylab = "Type I Error Rate",
             filename = "one_sample_typeI_error.pdf", hline = 0.05)

create_plots(list(power_t.test, power_wilcox.test, power_bootstrap.test),
             y_var = "power", title = "Power", ylab = "Power",
             filename = "one_sample_power.pdf", hline = 0.8)

# ----------------------------
# Two-sample case simulation
# ----------------------------

# Bootstrap two-sample test
bootstrap_two_sample_test <- function(x, y, B = 1000) {
  T_obs <- mean(x) - mean(y)
  combined <- c(x, y)
  n1 <- length(x)
  n2 <- length(y)
  
  x_centered <- x - mean(x) + mean(combined)
  y_centered <- y - mean(y) + mean(combined)
  
  boot_stats <- replicate(B, {
    x_star <- sample(x_centered, n1, replace = TRUE)
    y_star <- sample(y_centered, n2, replace = TRUE)
    mean(x_star) - mean(y_star)
  })
  
  p_val <- mean(abs(boot_stats) >= abs(T_obs))
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
        pval_boot_H0[i] <- bootstrap_two_sample_test(x, y, B)
        
        # power
        pval_t_H1[i] <- t.test(x, y + effect_size)$p.value
        pval_wilcox_H1[i] <- wilcox.test(x, y + effect_size)$p.value
        pval_boot_H1[i] <- bootstrap_two_sample_test(x, y + effect_size, B)
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
create_plots(list(TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test),
             y_var = "error", title = "Type I Error Rates", ylab = "Type I Error Rate",
             filename = "two_sample_typeI_error.pdf", hline = 0.05)

create_plots(list(power_t.test, power_wilcox.test, power_bootstrap.test),
             y_var = "power", title = "Power", ylab = "Power",
             filename = "two_sample_power.pdf", hline = 0.8)