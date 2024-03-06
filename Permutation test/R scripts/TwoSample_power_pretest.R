setwd("/Users/benedictkongyir/Desktop/OSU/Research/Permutation Test/New Simulations/R scripts")
source("User_defined_functions.R")

N <- 1e4;  alpha <- 0.05
dist_sum <- c("Standard Normal", "Uniform", "t", "Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", 
              "LogNormal", "Contaminated")
testvec <- c("KS", "SW", "JB", "DAP", "AD", "SF")
nvec <- c(8, 10, 15, 20, 30, 50)

set.seed(33)

powervec <- numeric(length(nvec) * length(dist_sum) * length(testvec))
power_t_test <- power_permutation <- power_pretest <- array(powervec, dim = c(6, 10, 6), 
                                                dimnames = list(nvec, dist_sum, testvec))

calculate_test_statistic <- function(x) {
  return(mean(x) / sqrt(var(x) / length(x)))
}

system.time({
  for (i in 1:length(nvec)) {
    n <- nvec[i]
    print(n)
    for (m in 1:length(testvec)) {
      test_type <- testvec[m]
      for (j in 1:length(dist_sum)) {
        dist <- dist_sum[j]
        pval_test_x <- pval_test_y <- numeric(N)
        for (k in 1:N) {
          x <- generate_data(n, dist)  
          y <- generate_data(n, dist)
          pval_test_x[k] <- generate_tests(x, test_type)$p.value
          pval_test_y[k] <- generate_tests(y, test_type)$p.value
        }
        power_pretest[i, j, m] <- mean(pval_test_x < alpha/2 | pval_test_y <alpha/2)
      }
    }
  }
})

power_pretest
save.image(paste0("Power_pretes_TwoSamples",".RData"))
