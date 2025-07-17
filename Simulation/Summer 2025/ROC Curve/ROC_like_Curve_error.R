# Read in user-defined functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

# Set working directory
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/Summer 2025/ROC Curve")

#-----------------------------------------------------------
# Parameters
alpha_pretest <- seq(from = 0.001, to = 0.05, by = 0.001)
n <- 10
Nsim <- 1e7
distributions <- c("Normal", "LogNormal")
effect_size <- 0.0
perm <- 1e3

#-----------------------------------------------------------
# Generate p-values function
generate_pval <- function(n, N, dist, effect_size, B) {
  pval_t.test <- pval_u.test <- pval_perm.test <- numeric(N)
  p_sw_x <- p_sw_y <- numeric(N)
  
  for (i in 1:N) {
    x <- generate_data(n, dist)
    y <- generate_data(n, dist)
    
    # SW-test p-values
    p_sw_x[i] <- shapiro.test(x)$p.value
    p_sw_y[i] <- shapiro.test(y)$p.value
    
    # downstream test p-values
    pval_t.test[i] <- t.test(x, y + effect_size)$p.value
    pval_u.test[i] <- wilcox.test(x, y + effect_size)$p.value
    #pval_perm.test[i] <- two_sample_permutation_test(x, y + effect_size, B)
  }
  
  # Return all p-values
  return(list(
    p_sw_x = p_sw_x,
    p_sw_y = p_sw_y,
    pval_t.test = pval_t.test,
    pval_u.test = pval_u.test
    #pval_perm.test = pval_perm.test
  ))
}
#-----------------------------------------------------------
# Store power values
Type_one_error <- list()
results <- list()

for (dist in distributions) {
  results[[dist]] <- generate_pval(n = n, N = Nsim, dist = dist, effect_size = effect_size, B = perm)
  
  # Calculate power
  Type_one_error[[dist]] <- list(
    Type_one_error_t.test = mean(results[[dist]]$pval_t.test < 0.05),
    Type_one_error_wilcoxon.test = mean(results[[dist]]$pval_u.test < 0.05)
    #power_perm.test = mean(results[[dist]]$pval_perm.test < 0.05)
  )
  
  for (j in seq_along(alpha_pretest)) {
    alpha <- alpha_pretest[j]
    
    adaptive_pvals_wilcox <- ifelse(
      results[[dist]]$p_sw_x > alpha & results[[dist]]$p_sw_y > alpha,
      results[[dist]]$pval_t.test,
      results[[dist]]$pval_u.test
    )
    
    # adaptive_pvals_perm <- ifelse(
    #   results[[dist]]$p_sw_x > alpha & results[[dist]]$p_sw_y > alpha,
    #   results[[dist]]$pval_t.test,
    #   results[[dist]]$pval_perm.test
    # )
    
    Type_one_error[[dist]]$adaptive_wilcox[j] <- mean(adaptive_pvals_wilcox < 0.05)
    #Type_one_error[[dist]]$adaptive_perm[j] <- mean(adaptive_pvals_perm < 0.05)
    Type_one_error[[dist]]$pr_sw_vec[j] <- mean(results[[dist]]$p_sw_x <= alpha 
                                               | results[[dist]]$p_sw_y <= alpha)
  }
}

#-----------------------------------------------------------
# Save results
save(
  results,
  Type_one_error,
  n,
  Nsim,
  distributions,
  alpha_pretest,
  file = "ROC_like_curve_Type_I_error.RData"
)

#-------------------------------------------------------------------------------------------
# Evaluate Efficiency Gains/Losses
EPG  <- Type_one_error$LogNormal$adaptive_wilcox - Type_one_error$LogNormal$Type_one_error_t.test
EPL  <- Type_one_error$Normal$adaptive_wilcox - Type_one_error$Normal$Type_one_error_t.test 

# # Calculate Point estimates
# EPG_lognormal = Type_one_error$LogNormal$Type_one_error_wilcoxon.test - Type_one_error$LogNormal$Type_one_error_t.test
# EPL_normal = Type_one_error$Normal$Type_one_error_t.test - Type_one_error$Normal$Type_one_error_wilcoxon.test 

#----------------------Plot results-----------------------------
# Save the  plots
pdf("comparison_error_inflation.pdf", width = 14, height = 8)
par(mfrow = c(1, 2))
plot(alpha_pretest, EPL, type = "l", col = "blue", ylab = "EPL", xlab = expression(alpha), main = "Expected Inflation of Type I error (Normal)")
plot(alpha_pretest, EPG, type = "l", col = "red", ylab = "EPG", xlab = expression(alpha), main = "Expected Inflation of Type I error (LogNormal)")
dev.off()

# Save the ROC like curve plot
pdf("ROC_like_curve_error.pdf", width = 6, height = 6)
plot(EPL, EPG, type = "l", col = "blue",
     xlab = "Expected Inflation of Type I error (Normal)",
     ylab = "Expected Inflation of Type I error (LogNormal)",
     main = "ROC like Curve: Inflation of Type I error")
#points(x = EPL_normal, y = EPG_lognormal)
dev.off()


