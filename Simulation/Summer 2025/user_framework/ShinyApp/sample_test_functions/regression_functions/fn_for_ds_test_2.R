fn_for_ds_test_2 <- function(data, n_boot = 1000) {
  if (!require(boot)) install.packages("boot")
  library(boot)
  
  # 1. Fit H0: y = μ + error  → get residuals from intercept‐only model
  mu_hat   <- mean(data$y)
  resid0   <- data$y - mu_hat
  df_boot  <- data.frame(x = data$x, resid0 = resid0)
  
  # 2. Observed slope
  observed_slope <- coef(lm(y ~ x, data = data))[2]
  
  # 3. Statistic: sample residuals with replacement, add to μ, refit slope
  boot_stat <- function(df, idx) {
    y_star <- mu_hat + df$resid0[idx]    # bootstrap y under H0
    model  <- lm(y_star ~ df$x)         # refit slope
    coef(model)[2]
  }
  
  # 4. Run bootstrap
  boot_results <- boot::boot(df_boot, statistic = boot_stat, R = n_boot)
  
  # 5. Two‐sided p‐value
  p_value <- mean(abs(boot_results$t) >= abs(observed_slope))
  
  return(list(p.value = p_value))
}
