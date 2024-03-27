source("/home/kongyir/spring2024/User_defined_functions.R")
N <- 1e4; P <- 1e4; alpha <- 0.05
dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Laplace","Chi-Square", "Gamma", "Weibull", 
              "LogNormal", "Contaminated")
nvec <- c(8, 10, 15, 20, 30, 50)
d.vec <- c(0.25, 0.5, 0.75)

set.seed(33)

powervec <- numeric(length(nvec) * length(dist_sum) * length(d.vec))
power_t_test <- power_permutation <- array(powervec, dim = c(6, 10, 3), 
                                           dimnames = list(nvec, dist_sum, d.vec))

# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}
system.time({
  for (i in 1:length(nvec)) {
    n <- nvec[i]
    print(n)
    for (m in 1:length(d.vec)) {
      d <- d.vec[m]
      for (j in 1:length(dist_sum)) {
        dist <- dist_sum[j]
        pval <- pval_perm <- numeric(N)
        for (k in 1:N) {
          x <- generate_data(n, dist) 
          y <- generate_data(n, dist) + d
          permuted_data <- c(x, y)
          pval[k] <- t.test(x, y, var.equal = F)$p.value
          
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
    }
  }
  powerloss <- power_permutation - power_t_test
})
power_t_test
power_permutation
powerloss

save.image(paste0("TwoSamples_powerloss",".RData"))
