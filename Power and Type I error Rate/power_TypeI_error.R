setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Power and Type I error Rate")
source("Defined functions.R")

N <- 1e3; alpha <- 0.05
nvec <- c(5, 8, 10, 20, 50)
d.vec <- c(0.25, 0.5, 0.75)
dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Chi-Square", "Gamma", "Weibull",  "LogNormal", "Contaminated")
powervec <- numeric(length(nvec) * length(dist_sum) * length(d.vec))
power_t_test <- array(powervec, dim = c(5, 9, 3), dimnames = list(nvec, dist_sum, d.vec))
set.seed(33)
system.time({
  for (i in 1:length(nvec)) {
    n <- nvec[i]
    print(n)
    for (m in 1:length(d.vec)) {
      d <- d.vec[m]
      for (j in 1:length(dist_sum)) {
        dist <- dist_sum[j]
        pval <- numeric(N).          
        for (k in 1:N) {
          x <- generate_data(n, dist) 
          y <- generate_data(n, dist) + d
          pval[k] <- t.test(x, y, var.equal = F)$p.value
        }
        power_t_test[i, j, m] <- mean(pval < alpha)
      }
    }
  }
})

power_t_test

save.image(paste0("TwoSamples_powerloss",".RData"))


# Load necessary libraries and source functions
source("OneSample functions.R")

# Set parameters
N <- 1e4
alpha <- 0.05
dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Chi-Square", "Gamma", "Weibull", 
              "LogNormal", "Contaminated")
testvec <- c("KS", "SW", "JB", "DAP")
nvec <- c(8, 10, 15, 20, 30, 50, 80)

# Initialize an empty data frame to store results
results <- data.frame(N = numeric(), Distribution = character(), Test = character(), Power = numeric())

# Set seed
set.seed(33)

# Loop over parameters and compute power for each combination
for (i in 1:length(nvec)) {
  n <- nvec[i]
  for (m in 1:length(testvec)) {
    test_type <- testvec[m]
    for (j in 1:length(dist_sum)) {
      dist <- dist_sum[j]
      pval_test <- numeric(N)
      for (k in 1:N) {
        x <- generate_data(n, dist)  
        pval_test[k] <- generate_tests(x, test_type)$p.value
      }
      power <- mean(pval_test < alpha)
      # Store results in the data frame
      results <- rbind(results, data.frame(N = n, Distribution = dist, Test = test_type, Power = power))
    }
  }
}

# View the results
print(results)
