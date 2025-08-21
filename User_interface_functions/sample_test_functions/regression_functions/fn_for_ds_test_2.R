fn_for_ds_test_2 <- function(data, n_boot = 1000) {
  # if (!require(boot)) install.packages("boot")
  # library(boot)
  mu_hat   <- mean(data$y)
  resid0   <- data$y - mu_hat
  df_boot  <- data.frame(x = data$x, resid0 = resid0)
  
  # Observed slope
  observed_slope <- coef(lm(y ~ x, data = data))[2]
  
  boot_stat <- function(df, idx) {
    y_star <- mu_hat + df$resid0[idx]    
    model  <- lm(y_star ~ df$x)         
    coef(model)[2]
  }
  
  boot_results <- boot::boot(df_boot, statistic = boot_stat, R = n_boot)
  
  p_value <- mean(abs(boot_results$t) >= abs(observed_slope))
  
  return(list(p.value = p_value))
}
