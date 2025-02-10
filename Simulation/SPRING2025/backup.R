# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# bootstrap function for one-sample location test
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bootstrap_one_sample <- function(x, effect_size, alpha , n_bootstrap, sample_size ) {
  pval <- c()
  for(j in seq_along(sample_size)){
    n <- sample_size[j]
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
    pval[j] <- mean(test_statistics_alt < lower_critical_value | test_statistics_alt > upper_critical_value)
  }
  
  return(pval)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# bootstrap function for two-sample location test
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# bootstrap function for two-sample location test
bootstrap_two_sample <- function(x1, x2, effect_size, alpha , n_bootstrap, sample_size ) {
  pval <- c()
  for(j in seq_along(sample_size)){
    n1 <- n2 <- sample_size[j]
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
    pval[j] <- mean(test_statistics_alt < lower_critical_value | test_statistics_alt > upper_critical_value)
  }
  
  return(pval)
}
