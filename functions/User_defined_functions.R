# load necessary libraries
library("nortest")
library("dgof")
library("dplyr")
library(moments)
library(tseries)
library(LaplacesDemon)
library(VGAM)
#Generate data from different distribution but located similarly
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
#%%%%%%%%%%% Apply different Normality tests %%%%%%%%%%%%%%%
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

# %%%%%%%%%%%% calculate one-sample  test statistic %%%%%%%%%%%%% 
OneSample_test_statistic <- function(x) {
  return((mean(x) * sqrt(length(x))) / sd(x))
}

# %%%%%%%%%%%% calculate two-sample  test statistic %%%%%%%%%%%%% 
TwoSample_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)) / (max(nvec) - min(nvec))
}

# OneSample Type I error Test Methods
OneSample_test <- function(x, test){
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
#TwoSample Type I error Test Methods
TwoSample_test <- function(x, y, test){
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

# %%%%%%%%%%%%%%%%%%%% POWER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# OneSample Power Test Methods
OneSample_Power.test <- function(x, test){
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
      index <- sample(c(-1, 1), length(x + d), replace = TRUE)
      sample_data <- index * abs(x + d)
      permuted_stat[j] <- OneSample_test_statistic(sample_data)
    }
    pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }
  return(pval)
}
#TwoSample Power Test Methods
TwoSample_Power.test <- function(x, y, test){
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

