
# Set the working directory (keep this line as is)
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Power and Type I error Rate")
source("OneSample functions.R")
# Load necessary functions (assuming they are defined elsewhere)

# Parameters
N <- 1e3
alpha <- 0.05
dist_sum <- c("Standard Normal", "Uniform", "t","Chi-Square", "Gamma", "Weibull", 
              "LogNormal")
shapiro_alpha <- c(0.1, 0.05, 0.01, 0.005, 0.000)
sample_size <- c(10, 20, 30, 40, 50)
set.seed(33)

# Initialize arrays
powervec <- numeric(length(dist_sum) * length(sample_size) * length(shapiro_alpha))
TypeI_error_ttest <- array(powervec, dim = c(7, 5, 5), dimnames = list(dist_sum, sample_size, shapiro_alpha))

# Simulation loop
for (i in 1:length(dist_sum)) {
  dist <- dist_sum[i]
  print(dist)
  for (m in 1:length(sample_size)) {
    n <- sample_size[m]
    for (j in 1:length(shapiro_alpha)) {
      pretest_alpha <- shapiro_alpha[j]
      num_sample <- 0
      pval_error <- numeric(N)
      while (num_sample < N) {
        x <- generate_data(n, dist)  
        y <- generate_data(n, dist)  
        if (shapiro.test(x)$p.value > pretest_alpha & shapiro.test(y)$p.value > pretest_alpha) {
          num_sample <- num_sample + 1
          pval_error[num_sample] <- t.test(x, y, var.equal = FALSE)$p.value
        }
      }
      TypeI_error_ttest[i, j, m] <- mean(pval_error < alpha)
    }
  }
}

# Print results or further analysis as needed
