# Load functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# Set working directory
setwd("/Users/benedictkongyir/Library/Mobile Documents/com~apple~CloudDocs/PhD Thesis/user_framework/ROC Curves")
# Parameters
alpha_pretest <- seq(from = 0.009, to = 0.1, by = 0.005)
sig_level <- seq(from = 0.01, to = 1, by = 0.01)
sample_size <- 10
Nsim <- 1e6
distributions <- c("Normal", "LogNormal")
effect_size <- 0.5


#------------------------------------------------------------
# Core Simulation Functions
#------------------------------------------------------------
generate_pval <- function(sample_size, N, dist, effect_size) {
  pval_t.test_H0 <- pval_wilcox.test_H0 <- numeric(N)
  pval_t.test_H1 <- pval_wilcox.test_H1 <- numeric(N)
  p_sw_x <- p_sw_y <- numeric(N)
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = N, style = 3) 
  
  for (i in 1:N) {
    x <- generate_data(sample_size, dist)
    y <- generate_data(sample_size, dist)
    
    # Pretest
    p_sw_x[i] <- shapiro.test(x)$p.value
    p_sw_y[i] <- shapiro.test(y)$p.value
    
    # Type I error p-values(under H0)
    pval_t.test_H0[i] <- t.test(x, y)$p.value
    pval_wilcox.test_H0[i] <- wilcox.test(x, y)$p.value
    
    # Power p-values (under H1)
    pval_t.test_H1[i] <- t.test(x, y + effect_size)$p.value
    pval_wilcox.test_H1[i] <- wilcox.test(x, y + effect_size)$p.value
    
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
#------------------------------------------------------------
#                     Perform Analysis  
# -----------------------------------------------------------
perform_analysis <- function(sample_size, N = Nsim, dist = dist, effect_size = effect_size, alpha_pretest){
  # Run simulation 
  results <- list()
  power_results <- list()
  error_results <- list()
  prob_sw_test <- list()
  
  for (dist in distributions) {
    # Generate data
    results[[dist]] <- generate_pval(sample_size, N = Nsim, dist = dist, effect_size = effect_size)
    
    # power Fixed alpha = 0.05 
      power_results[[dist]] <- list(
      power_t.test = mean(results[[dist]]$pval_t.test_H1 < 0.05),
      power_wilcox.test = mean(results[[dist]]$pval_wilcox.test_H1 < 0.05)
    )
    # Type I error Fixed alpha = 0.05 
      error_results[[dist]] <- list(
      error_t.test = mean(results[[dist]]$pval_t.test_H0 < 0.05),
      error_wilcox.test = mean(results[[dist]]$pval_wilcox.test_H0 < 0.05)
    )
    
    # Initialize vectors for adaptive test
    power_results[[dist]]$adaptive_wilcox <- numeric(length(alpha_pretest))
    error_results[[dist]]$adaptive_wilcox <- numeric(length(alpha_pretest))
    prob_sw_test[[dist]]$pr_sw_vec <- numeric(length(alpha_pretest))
    
    # Evaluate adaptive test at each alpha_pretest level
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      
      decision_region <- results[[dist]]$p_sw_x > alpha & results[[dist]]$p_sw_y > alpha
      
      # Adaptive Type I error
      adaptive_error_pvals <- ifelse(decision_region,
                                     results[[dist]]$pval_t.test_H0,
                                     results[[dist]]$pval_wilcox.test_H0)
      error_results[[dist]]$adaptive_wilcox[j] <- mean(adaptive_error_pvals < 0.05)
      
      # Adaptive Power
      adaptive_power_pvals <- ifelse(decision_region,
                                     results[[dist]]$pval_t.test_H1,
                                     results[[dist]]$pval_wilcox.test_H1)
      power_results[[dist]]$adaptive_wilcox[j] <- mean(adaptive_power_pvals < 0.05)
      # probability of passing sw test
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

#------------------------------------------------------------
# Analysis Functions
#------------------------------------------------------------

run_simulation <- function() {
  sim_data <- perform_analysis(distributions, sample_size, Nsim, effect_size, alpha_pretest)
  
  # Unpack for saving
  results <- sim_data$results
  power_results <- sim_data$power_results
  error_results <- sim_data$error_results
  
  save(
    results,
    power_results,
    error_results,
    alpha_pretest,
    sample_size,
    Nsim,
    distributions,
    file = "ROC_like_curve_for_power_error.RData"
  )
  
  return(sim_data)
}
#------------------------------------------------------------
# Compute ROC-like curve quantities
compute_roc_metrics <- function(power_results, error_results){
  # Calculate EPL and EPG
  EPG  <- power_results$LogNormal$adaptive_wilcox - power_results$LogNormal$power_t.test
  EPL  <- power_results$Normal$power_t.test - power_results$Normal$adaptive_wilcox
  
  # Calculate EITE and EDTE
  Expected_Inflation_lognormal <- error_results$LogNormal$adaptive_wilcox - error_results$LogNormal$error_t.test
  Expected_Inflation_normal    <- error_results$Normal$adaptive_wilcox - error_results$Normal$error_t.test
  
  # Point estimates (benchmark comparison)
  EPG_lognormal <- power_results$LogNormal$power_wilcox.test - power_results$LogNormal$power_t.test
  EPL_normal <- power_results$Normal$power_t.test - power_results$Normal$power_wilcox.test
  
  # Point estimates (benchmark comparison)
  Inflation_lognormal <- error_results$LogNormal$error_wilcox.test - error_results$LogNormal$error_t.test
  Inflation_normal <- error_results$Normal$error_wilcox.test - error_results$Normal$error_t.test
  
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
#                 Plot Power/Error Trade off
# ------------------------------------------------------------
plot_power_error_tradeoff <- function(alpha_pretest = alpha_pretest, metrics, file_name = "comparison_power_error.pdf") {
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

# calculate ROC metrics
metrics <- compute_roc_metrics(sim_data$power_results,sim_data$error_results)
# Create the plots
plot_power_error_tradeoff(
  alpha_pretest = alpha_pretest,
  metrics = metrics,
  file_name = "comparison_power_error.pdf"
)

# -------------------------------------------------------------------
#             Plot ROC-like Curves (EPG vs. EPL and Inflation)
# -------------------------------------------------------------------
plot_roc_like_curves <- function(metrics, file_name = "ROC_like_curve_combine.pdf") {
  pdf(file_name, width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  plot(metrics$EPL, metrics$EPG, type = "l", col = "blue", lwd = 2,
       xlab = "Power Loss (Normal)",
       ylab = "Power Gain (LogNormal)",
       main = "ROC-like Curve: EPG vs. EPL")
  legend("bottomright", legend = paste("Maximum Gain =", round(metrics$EPG_lognormal - max(metrics$EPG), 3)), title = "Optimal Point")
  
  plot(metrics$Expected_Inflation_normal, metrics$Expected_Inflation_lognormal, type = "l", col = "red", lwd = 2,
       xlab = "Type I Error Inflation (Normal)",
       ylab = "Type I Error Inflation (LogNormal)",
       main = "ROC-like Curve: Type I Inflation")
  legend("bottomright", legend = paste("Maximum Gain =", round(metrics$Inflation_lognormal - max(metrics$Expected_Inflation_lognormal), 3)), title = "Optimal Point")
  dev.off()
}

# create ROC like curves

plot_roc_like_curves(
  metrics = metrics,
  file_name = "ROC_like_curve_combine.pdf"
)

# ---------------------------------------------------------
# -------- Power vs Type I Error (ROC-like Curve) ---------
# ---------------------------------------------------------

power_Type_I_error_function <- function(sample_size, Nsim, distributions, effect_size, sig_levels) {
  roc_data <- data.frame()
  
  for (dist in distributions) {
    cat("Simulating for:", dist, "\n")
    res <- generate_pval(sample_size, N = Nsim, dist = dist, effect_size = effect_size)
    
    for (alpha in sig_levels) {
      
      # Compute power for each alpha
      power_t <- mean(res$pval_t.test_H1 < alpha)
      power_w <- mean(res$pval_wilcox.test_H1 < alpha)
      # compute Type I error for each alpha
      error_t <- mean(res$pval_t.test_H0 < alpha)
      error_w <- mean(res$pval_wilcox.test_H0 < alpha)
      
      # Append for t-test
      roc_data <- rbind(roc_data, data.frame(
        Distribution = dist,
        Method = "t-test",
        Alpha = alpha,
        Power = power_t,
        TypeIError = error_t
      ))
      
      # Append for Wilcoxon
      roc_data <- rbind(roc_data, data.frame(
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

# Run the simulation function
roc_data <- power_Type_I_error_function(
  sample_size = sample_size,
  Nsim = Nsim,
  distributions = distributions,
  effect_size = effect_size,
  sig_levels = sig_level
)

# Save results
save(
  sample_size,
  Nsim,
  distributions,
  effect_size,
  sig_level,
  roc_data,
  file = "Power_vs_TypeIError_ROC.RData"
)
# ---------------------------------------------------------------
#                Plot Power vs Type I Error ROC Curve 
# ---------------------------------------------------------------
pdf("Power_vs_TypeIError_ROC_By_Distribution.pdf", width = 10, height = 6)
plot_power_vs_error <- function(roc_data, file_name = "Power_vs_TypeIError_ROC.pdf") {
  pdf(file_name, width = 10, height = 6)
  par(mfrow = c(1, length(unique(roc_data$Distribution))))
  
  methods <- c("t-test", "Wilcoxon")
  colors <- c("blue", "red")
  pch_vals <- c(19, 17)
  
  for (dist in unique(roc_data$Distribution)) {
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Type I Error", ylab = "Power",
         main = paste("ROC-like Curve(", dist, ")", sep = ""))
    
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

# Run power vs Type I error
roc_data <- power_Type_I_error_function(
  sample_size = sample_size,
  Nsim = Nsim,
  distributions = distributions,
  effect_size = effect_size,
  sig_levels = sig_level
)

plot_power_vs_error(
  roc_data = roc_data,
  file_name = "Power_vs_TypeIError_ROC_curve.pdf"
)

