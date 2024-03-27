# Load necessary libraries
library(foreach)
library(doParallel)
# Source the file containing the required functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Permutation test/R scripts")
source("User_defined_functions.R")

# Specify the number of cores to use
num_cores <- 2  # Change this according to your CPU specs
# Register parallel backend
registerDoParallel(cores = num_cores)



N <- 1e3
P <- 1e3
alpha <- 0.05
dist_sum <- c("Standard Normal", "Uniform", "t")
nvec <- c(5, 8, 10, 20, 50)
d.vec <- c(0.25, 0.5, 0.75)

set.seed(33)

powervec <- numeric(length(nvec) * length(dist_sum) * length(d.vec))
power_t_test <- power_permutation <- array(powervec, dim = c(5, 3, 3), 
                                           dimnames = list(nvec, dist_sum, d.vec))

# Perform parallel processing
system.time({
  foreach(i = 1:length(nvec), .combine = 'c') %:%
    foreach(m = 1:length(d.vec), .combine = 'c') %:%
    foreach(j = 1:length(dist_sum), .combine = 'c') %dopar% {
      n <- nvec[i]
      d <- d.vec[m]
      dist <- dist_sum[j]
      print(n)
      pval <- pval_perm <- numeric(N)
      for (k in 1:N) {
        x <- generate_data(n, dist) 
        y <- generate_data(n, dist) + d
        permuted_data <- c(x, y)
        pval[k] <- t.test(x, y, var.equal = FALSE)$p.value
        
        observed_statistic <- calculate_test_statistic(x, y)
        
        permuted_statistics <- rep(0, P)
        for (l in 1:P) {
          sampled_data <- sample(permuted_data)
          permuted_data1 <- sampled_data[1:length(x)]
          permuted_data2 <- sampled_data[(length(x) + 1):(length(x) + length(y))]
          permuted_statistics[l] <- calculate_test_statistic(permuted_data1, permuted_data2)
        }
        pval_perm[k] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
      }
      power_t_test[i, j, m] <- mean(pval < alpha)
      power_permutation[i, j, m] <- mean(pval_perm < alpha)
    }
  powerloss <- power_permutation - power_t_test
})

# Combine results
power_t_test
power_permutation
powerloss

# Stop the parallel backend
stopImplicitCluster()
