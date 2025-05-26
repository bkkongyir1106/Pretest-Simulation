# Load required packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, ggplot2, dplyr, purrr, tidyr, e1071, tseries, nortest)

# Source your custom functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
# Define common parameters
sample_sizes_power <- c(8, 10, 15, 20, 25, 30, 50)
test_names <- c("t", "Wilcox", "t_Wilcox", "perm")

# Define Calculate_power in global scope
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

# bootstrap method
bootstrap_calulate <- function(twosamples = FALSE){
  if(twosamples){
    # Two sample case
    results <- bootstrap_two_sample(x1 = data_load$x, x2 = data_load$y, effect_size = effect_size, alpha = 0.05,  n_bootstrap =1e4, sample_size = sample_size)
  }else{
    # One sample case
    results <- bootstrap_two_sample(x = data_load$x, effect_size = effect_size, alpha = 0.05,  n_bootstrap =1e4, sample_size = sample_size)
  }
  return(results)
}

