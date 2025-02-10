bootstrap_one_sample_power <- function(x, effect_size = 0.5, alpha = 0.05, n_bootstrap =1e4 ) {
  n <- length(x)
  #  Compute critical values under the null hypothesis
  test_statistics_null <- numeric(n_bootstrap)  
  for (i in 1:n_bootstrap) {
    bootstrap_sample <- sample(x, n, replace = TRUE)  
    test_statistics_null[i] <- OneSample_test_statistic(bootstrap_sample)
  }
  lower_critical_value <- quantile(test_statistics_null, probs = alpha / 2)
  upper_critical_value <- quantile(test_statistics_null, probs = 1 - alpha / 2)
  
  # Simulate test statistics under the alternative hypothesis
  test_statistics_alt <- numeric(n_bootstrap)  
  for (i in 1:n_bootstrap) {
    bootstrap_sample <- sample(x, n, replace = TRUE)  
    bootstrap_sample_alt <- bootstrap_sample + effect_size  
    test_statistics_alt[i] <- OneSample_test_statistic(bootstrap_sample_alt)
  }
  
  # Compute power
  power <- mean(test_statistics_alt < lower_critical_value | test_statistics_alt > upper_critical_value)
  
  # Plot the results
  plot_data <- data.frame(
    Test_Statistic = c(test_statistics_null, test_statistics_alt),
    Distribution = factor(rep(c("Null", "Alternative"), each = n_bootstrap))
  )
  ggplot_object <- ggplot(plot_data, aes(x = Test_Statistic, fill = Distribution)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = lower_critical_value, linetype = "dashed", color = "red") +
    geom_vline(xintercept = upper_critical_value, linetype = "dashed", color = "red") +
    annotate("text", x = lower_critical_value, y = 0, label = "Lower Critical Value", vjust = -1, hjust = 1, color = "red") +
    annotate("text", x = upper_critical_value, y = 0, label = "Upper Critical Value", vjust = -1, hjust = 0, color = "red") +
    labs(
      title = "Bootstrap Distributions of Test Statistic",
      subtitle = paste("Power =", round(power, 3)),
      x = "Test Statistic",
      y = "Density"
    ) +
    theme_minimal()
  print(ggplot_object)
  return(power)
}
# Example usage:
set.seed(123) 
power_value <- bootstrap_one_sample_power(x = rexp(n = 100, rate = 1), effect_size = 0.5)
print(power_value)


# bootstrap function for two-sample location test
bootstrap_two_sample_power <- function(x1, x2, effect_size = 0.5, alpha = 0.05, n_bootstrap) {
  n1 <- length(x1)
  n2 <- length(x2)

  # Compute critical values under the null hypothesis
  test_statistics_null <- numeric(n_bootstrap)  
  for (i in 1:n_bootstrap) {
    # Resample under the null hypothesis 
    combined_sample <- c(x1, x2)
    bootstrap_sample1 <- sample(combined_sample, n1, replace = TRUE)
    bootstrap_sample2 <- sample(combined_sample, n2, replace = TRUE)
    test_statistics_null[i] <- TwoSample_test_statistic(bootstrap_sample1, bootstrap_sample2)
  }
  # Calculate critical values from the null distribution
  lower_critical_value <- quantile(test_statistics_null, probs = alpha / 2)
  upper_critical_value <- quantile(test_statistics_null, probs = 1 - alpha / 2)
  
  # Simulate test statistics under the alternative hypothesis
  test_statistics_alt <- numeric(n_bootstrap)  
  for (i in 1:n_bootstrap) {
    # Resample under the alternative hypothesis 
    bootstrap_sample1 <- sample(x1, n1, replace = TRUE)
    bootstrap_sample2 <- sample(x2, n2, replace = TRUE) + effect_size 
    test_statistics_alt[i] <- TwoSample_test_statistic(bootstrap_sample1, bootstrap_sample2)
  }
  
  # Compute power
  power <- mean(test_statistics_alt < lower_critical_value | test_statistics_alt > upper_critical_value)
  
  # Plot the results
  plot_data <- data.frame(
    Test_Statistic = c(test_statistics_null, test_statistics_alt),
    Distribution = factor(rep(c("Null", "Alternative"), each = n_bootstrap))
  )
  
  ggplot_object <- ggplot(plot_data, aes(x = Test_Statistic, fill = Distribution)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = lower_critical_value, linetype = "dashed", color = "red") +
    geom_vline(xintercept = upper_critical_value, linetype = "dashed", color = "red") +
    annotate("text", x = lower_critical_value, y = 0, label = "Lower Critical Value", vjust = -1, hjust = 1, color = "red") +
    annotate("text", x = upper_critical_value, y = 0, label = "Upper Critical Value", vjust = -1, hjust = 0, color = "red") +
    labs(
      title = "Bootstrap Distributions of Test Statistic (Two-Sample)",
      subtitle = paste("Power =", round(power, 3)),
      x = "Test Statistic",
      y = "Density"
    ) +
    theme_minimal()

  # display the plot
  print(ggplot_object)
  
  return(power)
}

bootstrap_two_sample_power(rnorm(100), rnorm(100),effect_size = 0.5, alpha = 0.05, n_bootstrap = 1000 )

