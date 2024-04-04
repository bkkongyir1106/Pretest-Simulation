rm(list = ls())
# Setting seed for reproducibility
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Parallel processing/parallel_sim/Type I error")
source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/User_defined_functions.R")
source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/utility.R")

{
  N = 1e4; treshold = 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Laplace", 
                "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated")
  nvec <- c(8, 10, 15, 20, 30, 50) 
  sig_level <- c(0.01, 0.05, 0.075)
  #set.seed(33)
}

powervec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
TypeI_error_t.test <- array(powervec,dim = c(length(nvec), length(dist_sum), 
             length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))

for(i in 1 : length(nvec))
  {
  n = nvec[i]
  print(n)
  for(j in 1 : length(sig_level))
  {
    alpha = sig_level[j]
    for(l in 1 : length(dist_sum))
    {
      dist = dist_sum[l]
      pval <- pval_perm <- numeric(N)
      for (k in 1:N) {
        x <- generate_data(n, dist) 
        y <- generate_data(n, dist) 
        pval[k] <- t.test(x, y)$p.value
      }
      TypeI_error_t.test[i, l, j] <- mean(pval < alpha)
    }
  }
}
TypeI_error_t.test