rm(list = ls())
# Install required packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest,gbm, tidyr, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Set working directory
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/User Interface")
# Call in external functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

sample_size = c(8, 10, 15, 20, 25, 30, 50 )
effect_size = 0.5
twosamples = FALSE
test = "t_Wilcox"
# --------------------------------------------------------------------------
# Perform Downstream Test( t, Wilcox, t/Wilcox, permutation)
# --------------------------------------------------------------------------
Calculate_power <- function(alpha, N, twosamples = FALSE, dist, sample_size, test, effect_size, B) {
  results <- numeric(length(sample_size))
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval <- numeric(N)
    for (i in 1:N) {
      if (twosamples) {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = TRUE)
        pval[i] <- TwoSample.test(data$x, data$y + effect_size, test, alpha, B)
      } else {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = FALSE)
        pval[i] <- OneSample.test(data$x + effect_size, test, alpha, B)
      }
    }
    results[j] <- mean(pval < alpha)
  }
  return(results)
}
# Produce Results
set.seed(12345)
results_test <- Calculate_power(alpha = 0.05, N = 1e4, twosamples = twosamples,
        dist = "Uniform", sample_size = sample_size, test = test, effect_size = effect_size, B = 1e3)
print(results_test)

# ------------------------------------------------------------------------
# Bootstrap Method
# ------------------------------------------------------------------------
# load data
data_load = Generate_data(datagen.type = 2, n = NULL, dist = NULL, two_samples = TRUE)

bootstrap_calulate <- function(twosamples = FALSE){
  if(twosamples){
    # Two sample case
    results <- bootstrap_two_sample(x1 = data_load$x, x2 = data_load$y, effect_size = effect_size, alpha = 0.05,  n_bootstrap =1e4, sample_size = sample_size)
    results2 <- bootstrap_two_sample_test(x = data_load$x, y = data_load$y, test = test, effect_size = effect_size, alpha= 0.05, n_bootstrap = 1e4, sample_size)
  }else{
    # One sample case
    results <- bootstrap_one_sample(x = data_load$x, effect_size = effect_size, alpha = 0.05,  n_bootstrap =1e4, sample_size = sample_size)
    results2 <- bootstrap_one_sample_test(x = data_load$x, test = test, effect_size = effect_size, alpha= 0.05, n_bootstrap = 1e4, sample_size)
  }
  return(list(results = results, results2 = results2))
}
# Produce Results
set.seed(12345)
Boot_results <- bootstrap_calulate(twosamples = twosamples)
print(Boot_results)

