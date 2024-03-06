setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Permutation functions")
source("OneSample functions.R")

N <- 1e4;  alpha <- 0.05
dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Chi-Square", "Gamma", "Weibull", 
              "LogNormal", "Contaminated")
testvec <- c("KS", "SW", "JB", "DAP")
nvec <- c(8, 10, 15, 20, 30, 50, 80)
d.vec <- c(0.01, 0.5, 0.75)

set.seed(33)

powervec <- numeric(length(nvec) * length(dist_sum) * length(testvec))
power_t_test <- power_permutation <- power_pretest <- array(powervec, dim = c(7, 9, 4), 
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
        pval_test<- numeric(N)
        for (k in 1:N) {
          x <- generate_data(n, dist)  
          pval_test[k] <- generate_tests(x, test_type)$p.value
        }
        power_pretest[i, j, m] <- mean(pval_test < alpha)
      }
    }
  }
})

power_pretest

save.image(paste0("OneSamples_pretest_power",".RData"))