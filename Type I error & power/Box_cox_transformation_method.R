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

# Box-Cox transformation function
box_cox_transform <- function(x) {
  # Find optimal lambda using MLE
  bc_result <- MASS::boxcox(x ~ 1, plotit = FALSE)
  lambda <- bc_result$x[which.max(bc_result$y)]
  
  # Apply transformation
  if (abs(lambda) < 1e-6) {
    transformed_x <- log(x)
  } else {
    transformed_x <- (x^lambda - 1) / lambda
  }
  
  return(transformed_x)
}

# Modified t-test with Box-Cox transformation if needed
t_test_with_boxcox <- function(x, mu0 = 0) {
  # Test for normality using Shapiro-Wilk
  shapiro_test <- shapiro.test(x)
  
  # If data is not normal, apply Box-Cox transformation
  if (shapiro_test$p.value < 0.05) {
    # Check if data is positive (required for Box-Cox)
    if (all(x > 0)) {
      x_transformed <- box_cox_transform(x)
      # Test transformed data for normality
      shapiro_transformed <- shapiro.test(x_transformed)
      
      # If transformation improves normality, use transformed data
      if (shapiro_transformed$p.value > shapiro_test$p.value) {
        x <- x_transformed
      }
    }
  }
  
  # Perform t-test
  t_test_result <- t.test(x, mu = mu0)
  return(t_test_result$p.value)
}

# Parameters
{
  Nsim <- 1e3
  alpha <- 0.05
  B <- 1e3
  effect_size <- 0.5  
  distributions <- c("Normal", "Uniform", "t", "Exponential") #, "Chi_Square", "LogNormal")
  sample_size <- c(8, 10, 15, 20, 25, 30, 50)
}

# Define plot function
create_plots <- function(results, y_var, title, ylab, filename, 
                         colors = c("blue", "green", "red"), 
                         pchs = c(16, 17, 21),
                         tests = c("t-test", "Wilcoxon test", "Bootstrap test"),
                         hline = NULL) {
  
  # pdf(filename, width = 10, height = 8)
  
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
  
  # dev.off()
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
                     .packages = c("LaplacesDemon", "VGAM", "MASS"),
                     .options.snow = opts) %:%
    foreach(dist = distributions) %dopar% {
      set.seed(12345)
      pval_t_H0 <- pval_wilcox_H0 <- pval_boot_H0 <- pval_t_boxcox_H0 <- numeric(Nsim)
      pval_t_H1 <- pval_wilcox_H1 <- pval_boot_H1 <- pval_t_boxcox_H1 <- numeric(Nsim)
      
      for (i in 1:Nsim) {
        x <- generate_data(n, dist)
        
        # Type I error
        pval_t_H0[i] <- t.test(x)$p.value
        pval_wilcox_H0[i] <- wilcox.test(x)$p.value
        pval_boot_H0[i] <- one_sample_bootstrap_test(x, mu0 = 0, B)
        pval_t_boxcox_H0[i] <- t_test_with_boxcox(x, mu0 = 0)
        
        # power
        pval_t_H1[i] <- t.test(x + effect_size)$p.value
        pval_wilcox_H1[i] <- wilcox.test(x + effect_size)$p.value
        pval_boot_H1[i] <- one_sample_bootstrap_test(x + effect_size, mu0 = 0, B)
        pval_t_boxcox_H1[i] <- t_test_with_boxcox(x + effect_size, mu0 = 0)
      }
      
      list(
        error_t_test = mean(pval_t_H0 < alpha),
        error_wilcox_test = mean(pval_wilcox_H0 < alpha),
        error_bootstrap_test = mean(pval_boot_H0 < alpha),
        error_t_boxcox_test = mean(pval_t_boxcox_H0 < alpha),
        power_t_test = mean(pval_t_H1 < alpha),
        power_wilcox_test = mean(pval_wilcox_H1 < alpha),
        power_bootstrap_test = mean(pval_boot_H1 < alpha),
        power_t_boxcox_test = mean(pval_t_boxcox_H1 < alpha)
      )
    }
  # close cluster
  close_cluster(my_cl)
  
  # Store results
  errorvec <- numeric(length(sample_size) * length(distributions))
  TypeI_error_t.test <- TypeI_error_wilcox.test <- TypeI_error_bootstrap.test <- TypeI_error_t_boxcox.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  power_t.test <- power_wilcox.test <- power_bootstrap.test <- power_t_boxcox.test <- 
    array(errorvec, dim = c(length(sample_size), length(distributions)),
          dimnames = list(sample_size, distributions))
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      TypeI_error_t.test[i, j] <- sim_out[[i]][[j]]$error_t_test
      TypeI_error_wilcox.test[i, j] <- sim_out[[i]][[j]]$error_wilcox_test
      TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
      TypeI_error_t_boxcox.test[i, j] <- sim_out[[i]][[j]]$error_t_boxcox_test
      power_t.test[i, j] <- sim_out[[i]][[j]]$power_t_test
      power_wilcox.test[i, j] <- sim_out[[i]][[j]]$power_wilcox_test
      power_bootstrap.test[i, j] <- sim_out[[i]][[j]]$power_bootstrap_test
      power_t_boxcox.test[i, j] <- sim_out[[i]][[j]]$power_t_boxcox_test
    }
  }
})

# Print results
cat("Type I error Rates for t-test\n"); print(TypeI_error_t.test)
cat("Type I error Rates for Wilcoxon-test\n"); print(TypeI_error_wilcox.test)
cat("Type I error Rates for Bootstrap test\n"); print(TypeI_error_bootstrap.test)
cat("Type I error Rates for t-test with Box-Cox\n"); print(TypeI_error_t_boxcox.test)
cat("Power for t-test\n"); print(power_t.test)
cat("Power for Wilcoxon-test\n"); print(power_wilcox.test)
cat("Power for Bootstrap test\n"); print(power_bootstrap.test)
cat("Power for t-test with Box-Cox\n"); print(power_t_boxcox.test)


# Create plots for one-sample case
# Create plots for Type I error rates
create_plots_typeI <- function() {
  # Layout: 2 rows x 2 columns for plots, 1 row for shared legend
  layout_matrix <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE)
  layout(mat = layout_matrix, heights = c(1, 1, 0.3))
  
  # Plotting setup
  par(mar = c(4, 4, 2, 1))
  colors <- c("blue", "green", "red", "purple")
  pchs <- c(16, 17, 21, 15)
  tests <- c("t-test", "Wilcoxon test", "Bootstrap test", "t-test with Box-Cox")
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get data for each test
    y_t <- TypeI_error_t.test[, j]
    y_w <- TypeI_error_wilcox.test[, j]
    y_b <- TypeI_error_bootstrap.test[, j]
    y_bc <- TypeI_error_t_boxcox.test[, j]
    
    # Determine ylim
    y_max <- max(y_t, y_w, y_b, y_bc, 0.2, na.rm = TRUE)
    
    # Create base plot
    plot(sample_size, y_t, 
         type = "b", 
         col = colors[1], 
         lty = 1, 
         pch = pchs[1],
         lwd = 2,
         ylim = c(0, y_max),
         xlab = "Sample Size", 
         ylab = "Type I Error Rate",
         main = dist_name)
    
    # Add other tests
    lines(sample_size, y_w, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_b, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
    lines(sample_size, y_bc, type = "b", col = colors[4], lty = 1, pch = pchs[4], lwd = 2)
    
    # Add reference line at alpha = 0.05
    abline(h = alpha, col = "gray", lty = 3, lwd = 2)
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
         horiz = FALSE,
         ncol = 2,
         bty = "n", 
         cex = 1.0)
}

# Create plots for Power
create_plots_power <- function() {
  # Layout: 2 rows x 2 columns for plots, 1 row for shared legend
  layout_matrix <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE)
  layout(mat = layout_matrix, heights = c(1, 1, 0.3))
  
  # Plotting setup
  par(mar = c(4, 4, 2, 1))
  colors <- c("blue", "green", "red", "purple")
  pchs <- c(16, 17, 21, 15)
  tests <- c("t-test", "Wilcoxon test", "Bootstrap test", "t-test with Box-Cox")
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Get data for each test
    y_t <- power_t.test[, j]
    y_w <- power_wilcox.test[, j]
    y_b <- power_bootstrap.test[, j]
    y_bc <- power_t_boxcox.test[, j]
    
    # Determine ylim
    y_max <- 1
    
    # Create base plot
    plot(sample_size, y_t, 
         type = "b", 
         col = colors[1], 
         lty = 1, 
         pch = pchs[1],
         lwd = 2,
         ylim = c(0, y_max),
         xlab = "Sample Size", 
         ylab = "Power",
         main = dist_name)
    
    # Add other tests
    lines(sample_size, y_w, type = "b", col = colors[2], lty = 1, pch = pchs[2], lwd = 2)
    lines(sample_size, y_b, type = "b", col = colors[3], lty = 1, pch = pchs[3], lwd = 2)
    lines(sample_size, y_bc, type = "b", col = colors[4], lty = 1, pch = pchs[4], lwd = 2)
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
         horiz = FALSE,
         ncol = 2,
         bty = "n", 
         cex = 1.0)
}

# Create PDF files with the plots
pdf("TypeI_error_rates_with_boxcox.pdf", width = 10, height = 8)
create_plots_typeI()
dev.off()

pdf("Power_rates_with_boxcox.pdf", width = 10, height = 8)
create_plots_power()
dev.off()

# Also display plots in R
par(mfrow = c(1, 1))  # Reset to single plot
create_plots_typeI()
create_plots_power()

# Create summary tables for easy comparison
cat("\n=== SUMMARY: AVERAGE TYPE I ERROR RATES ===\n")
avg_typeI <- data.frame(
  Distribution = distributions,
  t_test = colMeans(TypeI_error_t.test),
  Wilcoxon = colMeans(TypeI_error_wilcox.test),
  Bootstrap = colMeans(TypeI_error_bootstrap.test),
  BoxCox_t = colMeans(TypeI_error_t_boxcox.test)
)
print(avg_typeI)

cat("\n=== SUMMARY: AVERAGE POWER ===\n")
avg_power <- data.frame(
  Distribution = distributions,
  t_test = colMeans(power_t.test),
  Wilcoxon = colMeans(power_wilcox.test),
  Bootstrap = colMeans(power_bootstrap.test),
  BoxCox_t = colMeans(power_t_boxcox.test)
)
print(avg_power)

# Create individual comparison plots for each distribution
create_individual_plots <- function() {
  par(mfrow = c(2, 2))
  
  for (j in seq_along(distributions)) {
    dist_name <- distributions[j]
    
    # Type I error plot
    plot(sample_size, TypeI_error_t.test[, j], type = "b", col = "blue", lwd = 2,
         ylim = c(0, max(TypeI_error_t.test[, j], TypeI_error_wilcox.test[, j], 
                         TypeI_error_bootstrap.test[, j], TypeI_error_t_boxcox.test[, j], 0.2)),
         xlab = "Sample Size", ylab = "Type I Error Rate", 
         main = paste(dist_name, "- Type I Error"))
    lines(sample_size, TypeI_error_wilcox.test[, j], type = "b", col = "green", lwd = 2)
    lines(sample_size, TypeI_error_bootstrap.test[, j], type = "b", col = "red", lwd = 2)
    lines(sample_size, TypeI_error_t_boxcox.test[, j], type = "b", col = "purple", lwd = 2)
    abline(h = alpha, col = "gray", lty = 2)
    
    # Power plot
    plot(sample_size, power_t.test[, j], type = "b", col = "blue", lwd = 2,
         ylim = c(0, 1), xlab = "Sample Size", ylab = "Power", 
         main = paste(dist_name, "- Power"))
    lines(sample_size, power_wilcox.test[, j], type = "b", col = "green", lwd = 2)
    lines(sample_size, power_bootstrap.test[, j], type = "b", col = "red", lwd = 2)
    lines(sample_size, power_t_boxcox.test[, j], type = "b", col = "purple", lwd = 2)
  }
  
  # Add legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = c("t-test", "Wilcoxon", "Bootstrap", "Box-Cox t-test"),
         col = c("blue", "green", "red", "purple"), lwd = 2, lty = 1, 
         horiz = TRUE, bty = "n", cex = 0.8)
}

# Create individual comparison plots
pdf("Individual_comparisons_with_boxcox.pdf", width = 12, height = 10)
create_individual_plots()
dev.off()

# Display individual plots
create_individual_plots()
