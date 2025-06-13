source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

generate_pval <- function(n, N, dist, effect_size, B){
  pval_t.test <- pval_u.test <- pval_perm.test <- numeric(N)
  p_sw_x <- p_sw_y <- numeric(N)
  
  for(i in 1:N){
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

# Parameters
set.seed(12345)
alpha_pretest <- seq(from = 0.01, to = 0.1, by = 0.01)
sample_sizes <- c(8, 10, 20, 30, 40, 50)
Nsim <- 1e2
perm <- 1e2
distributions <- c("Normal", "LogNormal")
effect_size <- 0.5

# Store results
results_all <- list()

for (dist in distributions) {
  results_all[[dist]] <- list()
  for (n in sample_sizes) {
    results_all[[dist]][[n]] <- list()
    
    results <- generate_pval(n = n, N = Nsim, dist = dist, effect_size = effect_size, B = perm)

    power_t <- mean(results$pval_t.test < 0.05)
    power_u <- mean(results$pval_u.test < 0.05)
    power_perm <- mean(results$pval_perm.test < 0.05)
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      pval_adaptive <- ifelse(
        results$p_sw_x > alpha & results$p_sw_y > alpha,
        results$pval_t.test,
        results$pval_u.test
      )
      power_adaptive <- mean(pval_adaptive < 0.05)
      pr_sw = mean(results$p_sw_x <= alpha | results$p_sw_y <= alpha)
    }
    
    # Store average power over all alpha levels
    results_all[[dist]][[as.character(n)]] <- list(
      prob_sw = pr_sw,
      power_t_test = power_t,
      power_u_test = power_u,
      power_perm_test = power_perm,
      mean_power_adaptive = power_adaptive
    )
  }
}


# data frame per distribution
test_power_dfs <- list()

for (dist in distributions) {
  
  dist_results <- results_all[[dist]]
  
  probability_sw <- t_powers <- u_powers <- adaptive_powers <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    n <- as.character(sample_sizes[i])
    probability_sw[i] <- dist_results[[n]]$pr_sw
    t_powers[i] <- dist_results[[n]]$power_t_test
    u_powers[i] <- dist_results[[n]]$power_u_test
    adaptive_powers[i] <- dist_results[[n]]$mean_power_adaptive
  }
  
  # Create dataframe
  df <- data.frame(
    SampleSize = sample_sizes,
    prob_sw.test  = probability_sw,
    t_test = t_powers,
    wilcox_test = u_powers,
    adaptive_test = adaptive_powers
  )
  
  test_power_dfs[[dist]] <- df
}

# View result 
print(test_power_dfs$Normal)
print(test_power_dfs$LogNormal)

# plot power curves
df_normal <- test_power_dfs$Normal      
df_lognormal <- test_power_dfs$LogNormal

# Set up plots features
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))  

cols <- c("blue", "red", "darkgreen")
pchs <- c(16, 17, 15)
legend_labels <- c("t-test", "Wilcoxon", "Adaptive")
legend_title <- "Method"

# --------- Plot 1: Normal Distribution ----------
plot(df_normal$SampleSize, df_normal$t_test, type = "b", pch = pchs[1], col = cols[1],
     ylim = c(0, 1), xlab = "Sample Size", ylab = "Power", main = "Normal Distribution")

lines(df_normal$SampleSize, df_normal$wilcox_test, type = "b", pch = pchs[2], col = cols[2])
lines(df_normal$SampleSize, df_normal$adaptive_test, type = "b", pch = pchs[3], col = cols[3])

legend("bottomright", legend = legend_labels,title = legend_title , col = cols, pch = pchs, bty = "n")

# --------- Plot 2: LogNormal Distribution ----------
plot(df_lognormal$SampleSize, df_lognormal$t_test, type = "b", pch = pchs[1], col = cols[1],
     ylim = c(0, 1), xlab = "Sample Size", ylab = "Power", main = "LogNormal Distribution")

lines(df_lognormal$SampleSize, df_lognormal$wilcox_test, type = "b", pch = pchs[2], col = cols[2])
lines(df_lognormal$SampleSize, df_lognormal$adaptive_test, type = "b", pch = pchs[3], col = cols[3])

legend("bottomright", legend = legend_labels, title = legend_title , col = cols, pch = pchs, bty = "n")

# ROC curve plot

EPL <- df_normal$adaptive_test - df_normal$t_test
EPG <- df_lognormal$adaptive_test - df_lognormal$t_test

plot(EPL, EPG, type = "b", xlim = c(min(EPL), max(EPL)), ylim = c(min(EPG), max(EPG)),  col = "red", pch = pchs[3],
     xlab = "Normal", ylab = "Lognormal")


