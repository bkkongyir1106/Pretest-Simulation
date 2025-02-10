# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, dgof,  nortest, ggplot2, dplyr, LaplacesDemon, VGAM)

#Generate data from different distribution 
generate_data <- function(n, dist){
    if(dist == "Standard Normal"){ 
      x <- rnorm(n, mean = 0, sd = 1)
    }
    if(dist == "Chi-Square"){
      x <- (rchisq(n, df = 3) - 3)/sqrt(6) 
    }
    if(dist == "Gamma"){
      x <- (rgamma(n, shape = 3, rate = 0.1) - 30)/sqrt(300)
    }
    if(dist == "Exponential"){
      x <- rexp(n, rate = 1) - 1 
    }
    if(dist == "t"){
      x <- (rt(n, df = 7))/sqrt(7/5) 
    }
    if(dist == "Uniform"){
      x <- (runif(n, min = 0, max = 1) - 0.5)*sqrt(12) 
    }
    if(dist == "Laplace"){
      x <- rlaplace(n , location = 0, scale = 4)/sqrt(8)
    }
    if(dist == "Weibull"){
      x <- (rweibull(n, shape = 1, scale = 2) - 2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2))) 
    }
    if(dist == "LogNormal"){
      x <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1)) 
    }
   
    if(dist == "Contaminated"){
      br <- rbinom(n , size = 1 , prob = 0.75)
      sd_br <- sqrt(1 + br * 24)
      x <- rnorm(n, mean = 0, sd = sd_br)/sqrt(0.25 + 0.75*24)
    }
  if(dist == "Pareto"){
    shape = 3
    mean_pareto <- shape / (shape - 1)
    sd_pareto <- shape/(((shape - 1)^2)*(shape - 2))
    x <- (rpareto(n, shape = shape) - mean_pareto)/sd_pareto
  }
  return(x)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define the functions to generate data from either a specific, a mixture 
# distribution of load from a source file
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate_data <- function(datagen.type, n = NULL, dist = NULL, two_samples = FALSE, 
                          priors = NULL, x_weights = NULL, y_weights = NULL, ...) {
  
  # generation via function
  if (datagen.type == 1) {
    if (!is.null(priors) && length(priors) == 2 && !is.null(x_weights) && dist == "mixture") {
      x_weights <- x_weights / sum(x_weights)
      # Generate data for x 
      x <- vector(length = n)
      for (i in 1:n) {
        component <- sample(1:2, size = 1, prob = x_weights)
        x[i] <- generate_data(1, priors[component])
      }
      if (two_samples) {
        y_weights <- y_weights / sum(y_weights)
        # Generate data for y
        y <- vector(length = n)
        for (i in 1:n) {
          component <- sample(1:2, size = 1, prob = y_weights)
          y[i] <- generate_data(1, priors[component])
        }
      } else {
        y <- NULL  
      }
    } else {
      x <- generate_data(n, dist)  
      if (two_samples) {
        y <- generate_data(n, dist)  
      } else {
        y <- NULL  
      }
    }
  }
  # loading from a CSV file
  if (datagen.type == 2) {
    # Prompt user to select a CSV file
    file_path <- file.choose()  
    data <- read.csv(file_path, header = TRUE)  
    if (two_samples) {
      if (ncol(data) < 2) {
        stop("The CSV file should contain at least two columns for two samples (x and y).")
      }
      x <- data[[1]]  
      y <- data[[2]]  
    } else {
      x <- data[[1]]  
      y <- NULL       
    }
  }
  
  return(list(x = x, y = y))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define function to perform different Normality tests
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generate_tests <- function(x, test){
  if(test == "KS"){
    output <- lillie.test(x)
  }
  if(test == "SW"){
    output <- shapiro.test(x)
  }
  if(test == "JB"){
    output <- jarque.test(x)
  }
  if(test == "DAP"){
    output <- agostino.test(x)
  }
  if(test == "AD"){
    output <- ad.test(x)
  }
  if(test == "SF"){
    output <- sf.test(x)
  }
  if(test == "CVM"){
    output <- cvm.test(x)
  }
  if(test == "CHISQ"){
    output <- chisq.test(x)
  }
return(output)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define simple calculation functions 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# calculate one-sample  test statistic #
OneSample_test_statistic <- function(x) {
  return((mean(x) * sqrt(length(x))) / sd(x))
}
# calculate two-sample  test statistic #
TwoSample_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# compute the area under the curve 
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)) / (max(nvec) - min(nvec))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define functions for performing various downstream tests 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OneSample Type I error #
OneSample_TypeI_test <- function(x, test, alpha){
  if(test == "t"){
    pval <- t.test(x)$p.value
  }
  if(test == "Wilcox"){
    pval <- wilcox.test(x)$p.value
  }
  if(test == "t_Wilcox"){
    if(shapiro.test(x)$p.value > alpha){
      pval <- t.test(x)$p.value
    }else{
      pval <- wilcox.test(x)$p.value
    }
  }
  if(test == "perm"){
    # Perform permutation test
    observe_stat <- OneSample_test_statistic(x)
      permuted_stat <- numeric(B)
      for (j in 1:B) {
        index <- sample(c(-1, 1), length(x), replace = TRUE)
        sample_data <- index * abs(x)
        permuted_stat[j] <- OneSample_test_statistic(sample_data)
      }
      pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }
  return(pval)
}
# TwoSample Type I error #
TwoSample_TypeI_error_test <- function(x, y, test, alpha){
  if(test == "t"){
    pval <- t.test(x, y)$p.value
  }
  if(test == "Wilcox"){
    pval <- wilcox.test(x, y)$p.value
  }
  if(test == "t_Wilcox"){
    if(shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha){
      pval <- t.test(x, y)$p.value
    }else{
      pval <- wilcox.test(x, y)$p.value
    }
  }
  if(test == "perm"){
    # Perform permutation test
    data <- c(x, y)
    observe_stat <- TwoSample_test_statistic(x, y)
    permuted_stat <- numeric(B)
    for (j in 1:B) {
      sample_data <- sample(data)
      sample_x <- sample_data[1:length(x)]
      sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
      permuted_stat[j] <- TwoSample_test_statistic(sample_x, sample_y) 
    }
    pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }
  return(pval)
}

# Calculate OneSample Power #
d = 0.5
B = 1000
OneSample_Power.test <- function(x, test, alpha){
  if(test == "t"){
    pval <- t.test(x + d)$p.value
  }
  if(test == "Wilcox"){
    pval <- wilcox.test(x + d)$p.value
  }
  if(test == "t_Wilcox"){
    if(shapiro.test(x)$p.value > alpha){
      pval <- t.test(x + d)$p.value
    }else{
      pval <- wilcox.test(x + d)$p.value
    }
  }
  if(test == "perm"){
    # Perform permutation test
    observe_stat <- OneSample_test_statistic(x + d)
    permuted_stat <- numeric(B)
    for (j in 1:B) {
      index <- sample(c(-1, 1), length(x), replace = TRUE)
      sample_data <- index * abs(x + d)
      permuted_stat[j] <- OneSample_test_statistic(sample_data)
    }
    pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }
  return(pval)
}
# Calculate TwoSample Power #
TwoSample_Power.test <- function(x, y, test, alpha){
  if(test == "t"){
    pval <- t.test(x, y + d)$p.value
  }
  if(test == "Wilcox"){
    pval <- wilcox.test(x, y + d)$p.value
  }
  if(test == "t_Wilcox"){
    if(shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha){
      pval <- t.test(x, y + d)$p.value
    }else{
      pval <- wilcox.test(x, y + d)$p.value
    }
  }
  if(test == "perm"){
    # Perform permutation test
    data <- c(x, y + d)
    observe_stat <- TwoSample_test_statistic(x, y + d)
    permuted_stat <- numeric(B)
    for (j in 1:B) {
      sample_data <- sample(data)
      sample_x <- sample_data[1:length(x)]
      sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
      permuted_stat[j] <- TwoSample_test_statistic(sample_x, sample_y) 
    }
    pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }
  return(pval)
}

# bootstrap function for one-sample location test
bootstrap_one_sample_power <- function(x, effect_size, alpha , n_bootstrap, sample_size ) {
  power <- c()
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
    power[j] <- mean(test_statistics_alt < lower_critical_value | test_statistics_alt > upper_critical_value)
  }
  
  return(power)
}
# Example usage:
set.seed(123) 
power_value <- bootstrap_one_sample_power(x = rexp(n = 10, rate = 1), effect_size = 0.0, alpha = 0.05,  n_bootstrap =1e4, sample_size = c(10, 20, 30))
print(power_value)

# bootstrap function for two-sample location test
bootstrap_two_sample_power <- function(x1, x2, effect_size, alpha , n_bootstrap, sample_size ) {
  power <- c()
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
    power[j] <- mean(test_statistics_alt < lower_critical_value | test_statistics_alt > upper_critical_value)
  }
  
  return(power)
}

# Example usage:
set.seed(123) 
power_value <- bootstrap_two_sample_power(x1 = rexp(n = 10, rate = 1), x2 = rexp(n = 10, rate = 1), effect_size = 0.5, alpha = 0.05,  n_bootstrap =1e4, sample_size = c(10, 20, 30))
print(power_value)






