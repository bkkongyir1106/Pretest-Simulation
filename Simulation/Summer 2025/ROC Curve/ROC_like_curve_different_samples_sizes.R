source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# set wkdr:
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/Summer 2025/ROC Curve")
# ---------------------------------------------------------
generate_pval <- function(n, N, dist, effect_size, B) {
  pval_t.test <- pval_u.test <- pval_perm.test <- numeric(N)
  p_sw_x <- p_sw_y <- numeric(N)
  
  for (i in 1:N) {
    x <- generate_data(n, dist)
    y <- generate_data(n, dist)
    
    p_sw_x[i] <- shapiro.test(x)$p.value
    p_sw_y[i] <- shapiro.test(y)$p.value
    
    pval_t.test[i] <- t.test(x, y + effect_size)$p.value
    pval_u.test[i] <- wilcox.test(x, y + effect_size)$p.value
    pval_perm.test[i] <- two_sample_permutation_test(x, y + effect_size, B)
  }
  
  list(
    p_sw_x = p_sw_x,
    p_sw_y = p_sw_y,
    pval_t.test = pval_t.test,
    pval_u.test = pval_u.test,
    pval_perm.test = pval_perm.test
  )
}

# -----------------------------------------------------
# Parameters
set.seed(12345)
alpha_pretest <- seq(from = 0.01, to = 1, by = 0.01)
sample_sizes <- c(8, 10, 20, 30, 40, 50)
Nsim <- 1e4
perm <- 1e3
distributions <- c("Normal", "LogNormal")
effect_size <- 0.5
# -----------------------------------------------------
# Store results
results_all <- list()

for (dist in distributions) {
  results_all[[dist]] <- list()
  for (n in sample_sizes) {
    results_all[[dist]][[as.character(n)]] <- list()
    
    results <- generate_pval(n = n, N = Nsim, dist = dist, effect_size = effect_size, B = perm)
    
    power_t <- mean(results$pval_t.test < 0.05)
    power_u <- mean(results$pval_u.test < 0.05)
    power_perm <- mean(results$pval_perm.test < 0.05)
    
    # Store adaptive power for each alpha_pretest
    power_adaptive_vec <- numeric(length(alpha_pretest))
    pr_sw_vec <- numeric(length(alpha_pretest))
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      pval_adaptive <- ifelse(
        results$p_sw_x > alpha & results$p_sw_y > alpha,
        results$pval_t.test,
        results$pval_u.test
      )
      power_adaptive_vec[j] <- mean(pval_adaptive < 0.05)
      pr_sw_vec[j] <- mean(results$p_sw_x <= alpha | results$p_sw_y <= alpha)
    }
    
    results_all[[dist]][[as.character(n)]] <- list(
      power_t_test = power_t,
      power_u_test = power_u,
      power_perm_test = power_perm,
      adaptive_power = power_adaptive_vec,
      prob_reject_sw = pr_sw_vec
    )
  }
}

#-------------------------------------------------------
# save RData
save(
  results_all,
  sample_sizes,
  Nsim,
  perm,
  distributions,
  alpha_pretest,
  effect_size,
  file = "ROC_like_curve_different_n.RData"
)
#------------------------------------------------------

