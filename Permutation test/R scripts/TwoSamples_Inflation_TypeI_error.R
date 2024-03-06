setwd("/Users/benedictkongyir/Desktop/OSU/Research/Permutation Test/New Simulations/R scripts")
source("User_defined_functions.R")

N <- 1e3; P <- 1e4
dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Laplace","Chi-Square", "Gamma", "Weibull", 
              "LogNormal", "Contaminated")
nvec <- c(5, 8, 10, 20, 30, 50)
sig_level <- c(0.01, 0.05, 0.075)
set.seed(33)

powervec <- numeric(length(nvec) * length(sig_level) * length(dist_sum))
TypeI_error_t.test <- TypeI_error_perm.test <- array(powervec, dim = c(6, 10, 3), 
                                           dimnames = list(nvec, dist_sum, sig_level))

# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}
system.time({
  for (i in 1:length(nvec)) {
    n <- nvec[i]
    print(n)
    for (m in 1:length(sig_level)) {
      alpha <- sig_level[m]
      for (j in 1:length(dist_sum)) {
        dist <- dist_sum[j]
        pval <- pval_perm <- numeric(N)
        for (k in 1:N) {
          x <- generate_data(n, dist) 
          y <- generate_data(n, dist)
          data <- c(x, y)
          # t test
          pval[k] <- t.test(x, y, var.equal = F)$p.value
          # permutation
          observed_statistic <- calculate_test_statistic(x, y)
          permuted_statistics <- rep(0, P)
          for (l in 1:P) { 
            sampled_data <- sample(data)
            permuted_data1 <- sampled_data[1:length(x)]
            permuted_data2 <- sampled_data[(length(x) + 1):(length(x) + length(y))]
            permuted_statistics[l] <- calculate_test_statistic(permuted_data1, permuted_data2)
            
          }
          pval_perm[k] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
        }
        TypeI_error_t.test[i, j, m] <- mean(pval < alpha)
        TypeI_error_perm.test[i, j, m] <- mean(pval_perm < alpha)
      }
    }
  }
  inflation_error <- TypeI_error_t.test - TypeI_error_perm.test
})
TypeI_error_t.test
TypeI_error_perm.test
inflation_error
save.image(paste0("TwoSamples_Inflation_error",".RData"))
