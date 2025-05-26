# Load required packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, ggplot2, dplyr, purrr, tidyr, e1071, tseries, nortest)

# Source your custom functions
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/Rshiny/Version 5_20250308/user_defined_functions.R")
#source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
# Define common parameters
sample_sizes_power <- c(8, 10, 15, 20, 25, 30, 50)
test_names <- c("t", "Wilcox", "t_Wilcox", "perm")

# Define function to calculate power/Type I error in global scope
Calculate_power <- function(alpha, N, twosamples = FALSE, dist, sample_size, test, effect_size, B, custom_func_path = NULL) {
  results <- numeric(length(sample_size))
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval <- numeric(N)
    for (i in 1:N) {
      if (twosamples) {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = TRUE)
        if (test == "custom") {
          req(custom_func_path)
          source(custom_func_path, local = TRUE)
          pval[i] <- custom_two_sample(data$x, data$y)
        } else {
          pval[i] <- TwoSample.test(data$x, data$y , effect_size, test, alpha, B)
        }
      } else {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = FALSE)
        if (test == "custom") {
          req(custom_func_path)
          source(custom_func_path, local = TRUE)
          pval[i] <- custom_one_sample(data$x, effect_size)
        } else {
          pval[i] <- OneSample.test(data$x, effect_size, test, alpha, B)
        }
      }
    }
    results[j] <- mean(pval < alpha)
  }
  return(results)
}
