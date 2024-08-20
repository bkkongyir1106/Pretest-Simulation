# setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/COMBINE RESULTS_20240719/COMPARE METHODS/Area Under Curve")
# # set directories in local computer
# source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
# To run on cluster
setwd("/home/kongyir/spring2024/power")
source("/home/kongyir/spring2024/User_defined_functions.R")
source("/home/kongyir/spring2024/utility.R")

# loads the package 
#library("writexl") 
#library(openxlsx)

# Define parameters
set.seed(12345)
N <- 1e3 # Number of iterations
sample_sizes <- c(5, 10, 15, 20 , 25, 30, 50)
B <- 1e4 # Number of permutations
alpha <- 0.05  # Significance level
d <- 0.5 # effect size

distributions <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal") #, "Pareto")

# Function to calculate test statistic
calculate_test_statistic <- function(x) {
  return((mean(x)*sqrt(length(x))) / sqrt(var(x)))
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))/(max(sample_sizes) - min(sample_sizes))
}

# Define table to store results
powervec <- numeric(length(sample_sizes) * length(distributions))
power_t_table <- power_wilcox_table <- power_t_wilcox_table <- power_perm_table <- time_t_table <- time_wilcox_table <- 
  time_t_wilcox_table <- time_perm_table <- array(powervec,dim = c(length(sample_sizes), length(distributions)), dimnames = list(sample_sizes, distributions))

# set up progress bar
total_steps <- length(distributions) * length(sample_sizes) * N
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
progress <- 0

for (dist in distributions) {
  power_t_list <- power_wilcox_list <- power_t_wilcox_list <- power_perm_list <- power_shapiro_list <- numeric(length(sample_sizes))
  avg_time_t_list <- avg_time_wilcox_list <- avg_time_t_wilcox_list <- avg_time_perm_list <- numeric(length(sample_sizes))
  
  for (s in 1:length(sample_sizes)) {
    n <- sample_sizes[s]
    pval_t <- pval_wilcox <- pvals <- pval_perm <- numeric(N)
    time_t <- time_wilcox <- time_t_wilcox <- time_perm <- numeric(N)
    Shapiro_pval <- numeric(N)
    
    for (i in 1:N) {
      x <- generate_data(n, dist)

      # Perform t test
      time_t[i] <- system.time({
        pval_t[i] <- t.test(x + d)$p.value
      })
      
      # Wilcoxon test
      time_wilcox[i] <- system.time({
        pval_wilcox[i] <- wilcox.test(x + d)$p.value
      })
      
      # Perform t-test/Wilcoxon
      time_t_wilcox[i] <- system.time({
        if (shapiro.test(x)$p.value > alpha ) {
          Shapiro_pval[i] <- 1
          pvals[i] <- t.test(x + d)$p.value
        } else {
          pvals[i] <- wilcox.test(x + d)$p.value
        }
      })["elapsed"]
      
      # Perform permutation test
      observe_stat <- calculate_test_statistic(x + d)
      # Measure computation time for permutation test
      time_perm[i] <- system.time({
        permuted_stat <- numeric(B)
        for (j in 1:B) {
          index <- sample(c(-1, 1), length(x), replace = TRUE)
          sample_data <- index * abs(x + d)
          permuted_stat[j] <- calculate_test_statistic(sample_data) 
        }
        pval_perm[i] <- mean(abs(permuted_stat) >= abs(observe_stat))
      })["elapsed"]
    }
    
    # Update progress bar
    progress <- progress + 1
    setTxtProgressBar(pb, progress)
    
    # Calculate power
    power_t_list[s] <- mean(pval_t < alpha)
    power_wilcox_list[s] <- mean(pval_wilcox < alpha)
    power_t_wilcox_list[s] <- mean(pvals < alpha)
    power_perm_list[s] <- mean(pval_perm < alpha)
    
    # Calculate average computation time for each procedure
    avg_time_t_list[s] <- mean(time_t)
    avg_time_wilcox_list[s] <- mean(time_wilcox)
    avg_time_t_wilcox_list[s] <- mean(time_t_wilcox)
    avg_time_perm_list[s] <- mean(time_perm)
  }
  
  # Store results in tables
  power_t_table[, dist] <- power_t_list
  power_wilcox_table[, dist] <- power_wilcox_list
  power_t_wilcox_table[, dist] <- power_t_wilcox_list
  power_perm_table[, dist] <- power_perm_list
  time_t_table[, dist] <- avg_time_t_list
  time_wilcox_table[, dist] <- avg_time_wilcox_list
  time_t_wilcox_table[, dist] <- avg_time_t_wilcox_list
  time_perm_table[, dist] <- avg_time_perm_list
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))/(max(sample_sizes) - min(sample_sizes))
}

# Calculate areas
area_t <- apply(power_t_table, 2, compute_area, x = sample_sizes)
area_wilcox <- apply(power_wilcox_table, 2, compute_area, x = sample_sizes)
area_t_wilcox <- apply(power_t_wilcox_table, 2, compute_area, x = sample_sizes)
area_perm <- apply(power_perm_table, 2, compute_area, x = sample_sizes)

# Compute power differences
power_diff_t <- power_t_wilcox_table - power_t_table
power_diff_wilcox <- power_t_wilcox_table - power_wilcox_table
power_diff_perm <- power_t_wilcox_table - power_perm_table
# Compute time differences
time_diff_t <- time_t_wilcox_table - time_t_table
time_diff_wilcox <- time_t_wilcox_table - time_wilcox_table
time_diff_perm <- time_t_wilcox_table - time_perm_table
# close progress bar
close(pb)

# Print power tables
cat("Power for t test:\n")
print(power_t_table)
cat("\nPower for wilcoxon test:\n")
print(power_wilcox_table)
cat("Power for t-test/Wilcoxon test:\n")
print(power_t_wilcox_table)
cat("\nPower for permutation test:\n")
print(power_perm_table)

# Print computation time tables
cat("\nComputation time for t test:\n")
print(time_t_table)
cat("\nComputation time for wilcoxon test:\n")
print(time_wilcox_table)
cat("\nComputation time for t-test/Wilcoxon test:\n")
print(time_t_wilcox_table)
cat("\nComputation time for permutation test:\n")
print(time_perm_table)

# Print difference in power
cat("\nPower of two-stage - power of t test:\n")
print(power_diff_t)
cat("\nPower of two-stage - power of wilcoxon test:\n")
print(power_diff_wilcox)
cat("\nPower of two-stage - power of perm test:\n")
print(power_diff_perm)

# Print difference in time
cat("\nTime of two-stage - time of t test:\n")
print(time_diff_t)
cat("\nTime of two-stage - time of wilcoxon test:\n")
print(time_diff_wilcox)
cat("\nTime of two-stage - time of perm test:\n")
print(time_diff_perm)

# Print areas
cat("Area under the power curve for t-test:\n")
print(area_t)
cat("\nArea under the power curve for permutation test:\n")
print(area_wilcox)
cat("Area under the power curve for t-test/Wilcoxon:\n")
print(area_t_wilcox)
cat("\nArea under the power curve for permutation test:\n")
print(area_perm)

# Save data
save(sample_sizes, power_t_wilcox_table, power_perm_table, power_t_table, power_wilcox_table,
     time_t_table, time_wilcox_table, time_t_wilcox_table, time_perm_table,
     power_diff_t,power_diff_wilcox, power_diff_perm,
     time_diff_t, time_diff_wilcox, time_diff_perm,
     area_t, area_wilcox, area_t_wilcox, area_perm,
     file = "OneSampleAreaUnderCurve_power.RData" )
