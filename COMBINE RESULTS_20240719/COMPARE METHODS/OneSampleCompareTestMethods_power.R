setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Old stuff/COMPARISON OF DOWNSTREAM TEST METHODS")
# set directories in local computer
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
# To run on cluster
# setwd("/home/kongyir/spring2024/power")
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

# loads the package 
#library("writexl") 
#library(openxlsx)

#####################.Two-Stage Procedure vs. Permutation test########################
set.seed(12345)
N <- 1e1
sample_sizes <- c(5, 10, 15, 20, 25, 30, 50)
B <- 1e1 # Number of permutations
alpha <- 0.05  # Significance level
d <- 0.5
#distributions <- c("exponential", "lognormal", "chisquared", "pareto", "normal")
distributions <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal") #, "Pareto")


# Function to calculate test statistic
calculate_test_statistic <- function(x) {
  return((mean(x)*sqrt(length(x))) / sqrt(var(x)))
}


# Initialize tables to store results
power_t_table <- power_wilcox_table <- power_t_wilcox_table <- power_perm_table <- power_shapiro <- data.frame(matrix(NA, nrow = length(sample_sizes), ncol = length(distributions)))
time_t_table <- time_wilcox_table <- time_t_wilcox_table <- time_perm_table <- data.frame(matrix(NA, nrow = length(sample_sizes), ncol = length(distributions)))

# Set column names for the tables
colnames(power_t_table) <- colnames(power_wilcox_table) <- colnames(time_t_table) <- colnames(time_wilcox_table) <- colnames(power_shapiro) <- distributions
colnames(power_t_wilcox_table) <- colnames(power_perm_table) <- colnames(time_t_wilcox_table) <- colnames(time_perm_table) <- distributions

rownames(power_t_table) <- rownames(power_wilcox_table) <- rownames(time_t_table) <- rownames(time_wilcox_table) <- rownames(power_shapiro) <- sample_sizes
rownames(power_t_wilcox_table) <- rownames(power_perm_table) <- rownames(time_t_wilcox_table) <- rownames(time_perm_table) <- sample_sizes

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
      y <- generate_data(n, dist)
      
      # Perform t test
      time_t[i] <- system.time({
        pval_t[i] <- t.test(x, y + d)$p.value
      })
      
      # Wilcoxon test
      time_wilcox[i] <- system.time({
        pval_wilcox[i] <- wilcox.test(x, y + d)$p.value
      })
      
      # Perform t-test/Wilcoxon
      time_t_wilcox[i] <- system.time({
        if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
          Shapiro_pval[i] <- 1
          pvals[i] <- t.test(x, y + d)$p.value
        } else {
          pvals[i] <- wilcox.test(x, y + d)$p.value
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
    power_shapiro_list[s] <- 1-  mean(Shapiro_pval)
    
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
  power_shapiro[, dist] <- power_shapiro_list
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

# power of Shapiro Wilk test
cat("\nPower of Shapiro-Wilk test:\n")
print(1- power_shapiro)

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

# Print Normalized Computational Efficiency (NCE):
cat("\nNormalized Computational Efficiency (NCE): t test:\n")
print(time_t_table/time_t_table)
cat("\nNormalized Computational Efficiency (NCE): wilcoxon test:\n")
print(time_t_table/time_wilcox_table)
cat("\nNormalized Computational Efficiency (NCE): perm test:\n")
print(time_t_table/time_perm_table)

# Print areas
cat("Area under the power curve for t-test:\n")
print(area_t)
cat("\nArea under the power curve for permutation test:\n")
print(area_wilcox)
cat("Area under the power curve for t-test/Wilcoxon:\n")
print(area_t_wilcox)
cat("\nArea under the power curve for permutation test:\n")
print(area_perm)

close(pb)

save(sample_sizes, power_t_wilcox_table, power_perm_table, power_t_table, power_wilcox_table,
     time_t_table, time_wilcox_table, time_t_wilcox_table, time_perm_table,
     power_diff_t,power_diff_wilcox, power_diff_perm,
     time_diff_t, time_diff_wilcox, time_diff_perm,
     area_t, area_wilcox, area_t_wilcox, area_perm,
     power_shapiro,
     file = "OneSampleCompareTestMethods_power.RData" )
