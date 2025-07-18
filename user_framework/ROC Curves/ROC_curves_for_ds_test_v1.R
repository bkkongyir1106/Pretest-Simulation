# Load required packages
library(pbapply)

# Load user-defined functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

# Set working directory
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/user_framework/ROC Curves")

# Parameters
alpha_pretest <- seq(from = 0.009, to = 0.1, by = 0.005)
sig_level <- seq(from = 0.01, to = 1, by = 0.01)
sample_size <- 10
Nsim <- 1e2 
distributions <- c("Normal", "LogNormal")
effect_size <- 0.75
test_alpha <- 0.05  

# ------------------------------------------------------------
# Function to generate all required p-values 
generate_pval <- function(sample_size, N, dist, effect_size) {
  # Initialize vectors
  pval_t.test_H0 <- pval_wilcox.test_H0 <- numeric(N)
  pval_t.test_H1 <- pval_wilcox.test_H1 <- numeric(N)
  p_sw_x <- p_sw_y <- numeric(N)
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (i in 1:N) {
    x <- generate_data(sample_size, dist)
    y <- generate_data(sample_size, dist)
    
    # Pretest (Shapiro-Wilk)
    p_sw_x[i] <- shapiro.test(x)$p.value
    p_sw_y[i] <- shapiro.test(y)$p.value
    
    # Type I error p-values (under H0)
    pval_t.test_H0[i] <- t.test(x, y, var.equal = TRUE)$p.value
    pval_wilcox.test_H0[i] <- wilcox.test(x, y, exact = FALSE)$p.value
    
    # Power p-values (under H1)
    pval_t.test_H1[i] <- t.test(x, y + effect_size, var.equal = TRUE)$p.value
    pval_wilcox.test_H1[i] <- wilcox.test(x, y + effect_size, exact = FALSE)$p.value
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(
    p_sw_x = p_sw_x,
    p_sw_y = p_sw_y,
    pval_t.test_H0 = pval_t.test_H0,
    pval_wilcox.test_H0 = pval_wilcox.test_H0,
    pval_t.test_H1 = pval_t.test_H1,
    pval_wilcox.test_H1 = pval_wilcox.test_H1
  ))
}

# ------------------------------------------------------------
# Function to perform complete analysis
perform_analysis <- function(sample_size, N, distributions, effect_size, alpha_pretest) {
  results <- list()
  power_results <- list()
  error_results <- list()
  prob_sw_test <- list()
  
  for (dist in distributions) {
    cat("Processing distribution:", dist, "\n")
    results[[dist]] <- generate_pval(sample_size, N, dist, effect_size)
    
    # Calculate fixed test results
    power_results[[dist]] <- list(
      power_t.test = mean(results[[dist]]$pval_t.test_H1 < test_alpha),
      power_wilcox.test = mean(results[[dist]]$pval_wilcox.test_H1 < test_alpha)
    )
    
    error_results[[dist]] <- list(
      error_t.test = mean(results[[dist]]$pval_t.test_H0 < test_alpha),
      error_wilcox.test = mean(results[[dist]]$pval_wilcox.test_H0 < test_alpha)
    )
    
    # Initialize adaptive test results
    power_results[[dist]]$adaptive_wilcox <- numeric(length(alpha_pretest))
    error_results[[dist]]$adaptive_wilcox <- numeric(length(alpha_pretest))
    prob_sw_test[[dist]]$pr_sw_vec <- numeric(length(alpha_pretest))
    
    # Evaluate adaptive test at each alpha_pretest level
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      
      # Decision rule: use t-test if both samples appear normal
      decision_region <- results[[dist]]$p_sw_x > alpha & results[[dist]]$p_sw_y > alpha
      
      # Adaptive test p-values
      adaptive_error_pvals <- ifelse(decision_region,
                                     results[[dist]]$pval_t.test_H0,
                                     results[[dist]]$pval_wilcox.test_H0)
      
      adaptive_power_pvals <- ifelse(decision_region,
                                     results[[dist]]$pval_t.test_H1,
                                     results[[dist]]$pval_wilcox.test_H1)
      
      # Store results
      error_results[[dist]]$adaptive_wilcox[j] <- mean(adaptive_error_pvals < test_alpha)
      power_results[[dist]]$adaptive_wilcox[j] <- mean(adaptive_power_pvals < test_alpha)
      prob_sw_test[[dist]]$pr_sw_vec[j] <- mean(!decision_region)
    }
  }
  
  return(list(
    results = results,
    power_results = power_results,
    error_results = error_results,
    prob_sw_test = prob_sw_test
  ))
}

# ------------------------------------------------------------
# Function to generate ROC data for power vs type I error
generate_roc_data <- function(sample_size, N, distributions, effect_size, sig_levels) {
  
  roc_data <- data.frame()
  
  for (dist in distributions) {
    cat("Generating ROC data for:", dist, "\n")
    res <- generate_pval(sample_size, N, dist, effect_size)
    
    for (alpha in sig_levels) {
      # Compute power & Type I error
      power_t <- mean(res$pval_t.test_H1 < alpha)
      power_w <- mean(res$pval_wilcox.test_H1 < alpha)
      error_t <- mean(res$pval_t.test_H0 < alpha)
      error_w <- mean(res$pval_wilcox.test_H0 < alpha)
      
      # Append results
      roc_data <- rbind(roc_data, 
                        data.frame(
                          Distribution = dist,
                          Method = "t-test",
                          Alpha = alpha,
                          Power = power_t,
                          TypeIError = error_t
                        ),
                        data.frame(
                          Distribution = dist,
                          Method = "Wilcoxon",
                          Alpha = alpha,
                          Power = power_w,
                          TypeIError = error_w
                        ))
    }
  }
  
  return(roc_data)
}

# ------------------------------------------------------------
# Function to compute ROC-like metrics
compute_roc_metrics <- function(power_results, error_results) {
  EPG <- power_results$LogNormal$adaptive_wilcox - power_results$LogNormal$power_t.test
  EPL <- power_results$Normal$power_t.test - power_results$Normal$adaptive_wilcox
  
  Expected_Inflation_lognormal <- error_results$LogNormal$adaptive_wilcox - error_results$LogNormal$error_t.test
  Expected_Inflation_normal <- error_results$Normal$adaptive_wilcox - error_results$Normal$error_t.test
  
  # Point estimates for benchmark comparison
  EPG_lognormal <- power_results$LogNormal$power_wilcox.test - power_results$LogNormal$power_t.test
  EPL_normal <- power_results$Normal$power_t.test - power_results$Normal$power_wilcox.test
  
  return(list(
    EPL = EPL,
    EPG = EPG,
    Expected_Inflation_lognormal = Expected_Inflation_lognormal,
    Expected_Inflation_normal = Expected_Inflation_normal,
    EPG_lognormal = EPG_lognormal,
    EPL_normal = EPL_normal
  ))
}

# ------------------------------------------------------------
# Plotting Functions
# ------------------------------------------------------------

plot_power_error_tradeoff <- function(alpha_pretest, metrics, file_name) {
  pdf(file_name, width = 12, height = 10)
  par(mfrow = c(2, 2))
  
  plot(alpha_pretest, metrics$EPL, type = "l", col = "blue", lwd = 2,
       ylab = "Expected Power Loss (EPL)", xlab = expression(alpha),
       main = "Power Loss (Normal)")
  
  plot(alpha_pretest, metrics$EPG, type = "l", col = "red", lwd = 2,
       ylab = "Expected Power Gain (EPG)", xlab = expression(alpha),
       main = "Power Gain (LogNormal)")
  
  plot(alpha_pretest, metrics$Expected_Inflation_normal, type = "l", col = "orange", lwd = 2,
       ylab = "Type I Error Inflation", xlab = expression(alpha),
       main = "Inflation (Normal)")
  
  plot(alpha_pretest, metrics$Expected_Inflation_lognormal, type = "l", col = "green", lwd = 2,
       ylab = "Type I Error Inflation", xlab = expression(alpha),
       main = "Inflation (LogNormal)")
  
  dev.off()
}

plot_roc_like_curves <- function(metrics, file_name) {
  pdf(file_name, width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  plot(metrics$EPL, metrics$EPG, type = "l", col = "blue", lwd = 2,
       xlab = "Power Loss (Normal)", 
       ylab = "Power Gain (LogNormal)",
       main = "ROC-like Curve: EPG vs. EPL")
  legend("bottomright", legend = paste("EPG Gain =", round(metrics$EPG_lognormal, 3)), title = "Point Estimate")
  
  plot(metrics$Expected_Inflation_normal, metrics$Expected_Inflation_lognormal, type = "l", col = "red", lwd = 2,
       xlab = "Type I Error Inflation (Normal)",
       ylab = "Type I Error Inflation (LogNormal)",
       main = "ROC-like Curve: Type I Inflation")
  
  dev.off()
}

plot_power_vs_error_ROC_curve <- function(roc_data, file_name) {
  pdf(file_name, width = 10, height = 6)
  par(mfrow = c(1, length(unique(roc_data$Distribution))))
  
  methods <- c("t-test", "Wilcoxon")
  colors <- c("blue", "red")
  pch_vals <- c(19, 17)
  
  for (dist in unique(roc_data$Distribution)) {
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Type I Error", ylab = "Power",
         main = paste("ROC-like Curve (", dist, ")", sep = ""))
    
    for (m in seq_along(methods)) {
      method <- methods[m]
      data_subset <- subset(roc_data, Distribution == dist & Method == method)
      lines(data_subset$TypeIError, data_subset$Power, type = "l",
            col = colors[m], lwd = 2, pch = pch_vals[m])
    }
    
    legend("bottomright", legend = methods, col = colors, lwd = 2, pch = pch_vals, title = "Method")
  }
  
  dev.off()
}

# ------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------

# Run main simulation
sim_data <- perform_analysis(sample_size, Nsim, distributions, effect_size, alpha_pretest)

# Compute ROC metrics
metrics <- compute_roc_metrics(sim_data$power_results, sim_data$error_results)

# Generate ROC data
roc_data <- generate_roc_data(
  sample_size = sample_size,
  N = Nsim,
  distributions = distributions,
  effect_size = effect_size,
  sig_levels = sig_level
)

# Save all results
save(
  list = c("sim_data", "metrics", "roc_data", "alpha_pretest", "sample_size", 
           "Nsim", "distributions", "effect_size", "sig_level", "test_alpha"),
  file = "Complete_ROC_Analysis.RData"
)

# Generate all plots
plot_power_error_tradeoff(
  alpha_pretest = alpha_pretest,
  metrics = metrics,
  file_name = "comparison_power_error.pdf"
)

plot_roc_like_curves(
  metrics = metrics,
  file_name = "ROC_like_curve_combine.pdf"
)

plot_power_vs_error_ROC_curve(
  roc_data = roc_data,
  file_name = "Power_vs_TypeIError_ROC_By_Distribution.pdf"
)