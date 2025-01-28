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

## Normality tests Methods ###
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

# Bootstrap: onesample Type I error
n_bootstrap = 1e4
bootstrap_Onesample.TypeI_error.test <- function(x) {
  n1 <- length(x)
  test.statistics0 <- OneSample_test_statistic(x)
  test.statistics1 <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    sample1 <- sample(x, n1, replace = TRUE)
    test.statistics1[i] <- OneSample_test_statistic(sample1)
  }
  pval <- mean(abs(test.statistics1) >= abs(test.statistics0))
  return(pval)
}

# Bootstrap: one sample power
n_bootstrap = 1e4
d = 0.5
bootstrap_Onesample.power.test <- function(x) {
  n1 <- length(x)
  test.statistics0 <- OneSample_test_statistic(x + d)
  test.statistics1 <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    sample1 <- sample(x + d, n1, replace = TRUE)
    test.statistics1[i] <- OneSample_test_statistic(sample1)
  }
  
  pval <- mean(abs(test.statistics1) >= abs(test.statistics0))
  return(pval)
}

# Bootstrap: two samples Type I error
n_bootstrap = 1e4
bootstrap_twosample.TypeI_error.test <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  test.statistics0 <- TwoSample_test_statistic(x, y)
  test.statistics1 <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    sample1 <- sample(x, n1, replace = TRUE)
    sample2 <- sample(y, n2, replace = TRUE)
    test.statistics1[i] <- TwoSample_test_statistic(sample1, sample2)
  }
  
  pval <- mean(abs(test.statistics1) >= abs(test.statistics0))
  return(pval)
}

# Bootstrap: two samples Power
n_bootstrap = 1e4
d = 0.5
bootstrap_twosample.power.test <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  test.statistics0 <- TwoSample_test_statistic(x, y + d)
  test.statistics1 <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    sample1 <- sample(x, n1, replace = TRUE)
    sample2 <- sample(y + d, n2, replace = TRUE)
    test.statistics1[i] <- TwoSample_test_statistic(sample1, sample2)
  }
  
  pval <- mean(abs(test.statistics1) >= abs(test.statistics0))
  return(pval)
}
# example
bootstrap_twosample.TypeI_error.test(x = rnorm(100), y = rnorm(100))

# Bootstrap: two samples Type I error
n_bootstrap = 1e4
bootstrap_twosample.TypeI_error.test <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  test.statistics0 <- TwoSample_test_statistic(x, y)
  test.statistics1 <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    sample1 <- sample(x, n1, replace = TRUE)
    sample2 <- sample(y, n2, replace = TRUE)
    test.statistics1[i] <- TwoSample_test_statistic(sample1, sample2)
  }
  
  pval <- mean(abs(test.statistics1) >= abs(test.statistics0))
  return(pval)
}
