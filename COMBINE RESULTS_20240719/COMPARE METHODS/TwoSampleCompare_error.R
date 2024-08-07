# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research.05.10.2024/Power")
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

#####################.Two-Stage Procedure vs. Permutation test########################
set.seed(12345)
N <- 1e3
sample_sizes <- c(5, 10, 15, 20, 25, 30, 50)
B <- 1e4 # Number of permutations
alpha <- 0.05  # Significance level
d <- 0.5
#distributions <- c("exponential", "lognormal", "chisquared", "pareto", "normal")
distributions <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal") #, "Pareto")

# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))/max(sample_sizes)
}

# Initialize tables to store results
error_t_table <- error_wilcox_table <- error_t_wilcox_table <- error_perm_table  <- power_shapiro <- data.frame(matrix(NA, nrow = length(sample_sizes), ncol = length(distributions)))
time_t_table <- time_wilcox_table<- time_t_wilcox_table <- time_perm_table <- data.frame(matrix(NA, nrow = length(sample_sizes), ncol = length(distributions)))

# Setting column names for the tables
colnames(error_t_table) <- colnames(error_wilcox_table) <-colnames(time_t_table) <-colnames(time_wilcox_table) <-colnames(power_shapiro)<- distributions
colnames(error_t_wilcox_table) <- colnames(error_perm_table) <-colnames(time_t_wilcox_table) <-colnames(time_perm_table) <- distributions

rownames(error_t_table) <- rownames(error_wilcox_table) <- rownames(time_t_table) <- rownames(time_wilcox_table) <-rownames(power_shapiro)<- sample_sizes
rownames(error_t_wilcox_table) <- rownames(error_perm_table) <- rownames(time_t_wilcox_table) <- rownames(time_perm_table) <- sample_sizes

total_steps <- length(distributions) * length(sample_sizes) * N
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
progress <- 0

for (dist in distributions) {
  error_t_list <- error_wilcox_list <- error_t_wilcox_list <- error_perm_list <- power_shapiro_list<-  numeric(length(sample_sizes))
  avg_time_t_list <-avg_time_wilcox_list <- avg_time_t_wilcox_list <-avg_time_perm_list <-numeric(length(sample_sizes)) 
  
  for (s in 1:length(sample_sizes)) {
    n <- sample_sizes[s]
    # Initialize empty vectors to store values
    pval_t <- pval_wilcox <- pvals <- pval_perm<- Shapiro_pval <- numeric(N)
    time_t <- time_wilcox <- time_t_wilcox <- time_perm <- numeric(N)
    for (i in 1:N) {
      x <- generate_data(n, dist)
      y <- generate_data(n, dist) 
      
      # perform t test
      time_t[i] <- system.time(
        {
          pval_t[i] <- t.test(x, y)$p.value
        }
      )
      # Wilcoxon test
      time_wilcox[i] <- system.time(
        {
          pval_wilcox[i] <- wilcox.test(x, y)$p.value
        }
      )
      # MPerform t-test/Wilcoxon
      time_t_wilcox[i] <- system.time({
        if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
          pvals[i] <- t.test(x, y )$p.value
          Shapiro_pval[i] <- 1
        } else {
          pvals[i] <- wilcox.test(x, y)$p.value
        }
      })["elapsed"]
      
      data <- c(x, y)
      observe_stat <- calculate_test_statistic(x, y)
      
      # Perform permutation test
      time_perm[i] <- system.time({
        permuted_stat <- numeric(B)
        for (j in 1:B) {
          sample_data <- sample(data)
          sample_x <- sample_data[1:length(x)]
          sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
          permuted_stat[j] <- calculate_test_statistic(sample_x, sample_y) 
        }
        pval_perm[i] <- mean(abs(permuted_stat) >= abs(observe_stat))
      })["elapsed"]
    }
    
    # Calculate areas
    area_t <- apply(error_t_table - 0.05, 2, compute_area, x = sample_sizes)
    area_wilcox <- apply(error_wilcox_table - 0.05, 2, compute_area, x = sample_sizes)
    area_t_wilcox <- apply(error_t_wilcox_table - 0.05, 2, compute_area, x = sample_sizes)
    area_perm <- apply(error_perm_table - 0.05, 2, compute_area, x = sample_sizes)
    
    # Update progress bar
    progress <- progress + 1
    setTxtProgressBar(pb, progress)
    # power of Shapiro-Wilk test
    power_shapiro_list[s] <- 1 - mean(Shapiro_pval)
    
    # Calculate power 
    error_t_wilcox_list[s] <- mean(pvals < alpha)
    error_perm_list[s] <- mean(pval_perm < alpha)
    error_t_list[s] <- mean(pval_t < alpha)
    error_wilcox_list[s] <- mean(pval_wilcox < alpha)
    
    # Calculate average computation time for each procedure
    avg_time_t_wilcox_list[s] <- mean(time_t_wilcox)
    avg_time_perm_list[s] <- mean(time_perm)
    avg_time_t_list[s] <- mean(time_t)
    avg_time_wilcox_list[s] <- mean(time_wilcox)
  }
  
  # Store results in tables
  power_shapiro[, dist] <- power_shapiro_list
  error_t_table[, dist] <- error_t_list
  error_wilcox_table[, dist] <- error_wilcox_list
  error_t_wilcox_table[, dist] <- error_t_wilcox_list
  error_perm_table[, dist] <- error_perm_list
  time_t_table[, dist] <- avg_time_t_list
  time_wilcox_table[, dist] <- avg_time_wilcox_list
  time_t_wilcox_table[, dist] <- avg_time_t_wilcox_list
  time_perm_table[, dist] <- avg_time_perm_list
  
}

NCE_t <- time_t_table/time_t_table
NCE_wilcox <- time_t_table/time_wilcox_table
NCE_t_wilcox <-time_t_table/time_t_wilcox_table
NCE_perm <- time_t_table/time_perm_table
close(pb)

# Print Type I error
cat("Power for t test:\n")
print(error_t_table)
cat("\nPower for wilcoxon test:\n")
print(error_wilcox_table)
cat("Power for t-test/Wilcoxon test:\n")
print(error_t_wilcox_table)
cat("\nPower for permutation test:\n")
print(error_perm_table)

# Print computation time tables
cat("\nComputation time for t test:\n")
print(time_t_table)
cat("\nComputation time for wilcoxon test:\n")
print(time_wilcox_table)
cat("\nComputation time for t-test/Wilcoxon test:\n")
print(time_t_wilcox_table)
cat("\nComputation time for permutation test:\n")
print(time_perm_table)

# Print Normalized Computational Efficiency (NCE):
cat("\nNormalized Computational Efficiency (NCE): t test:\n")
print(NCE_t)
cat("\nNormalized Computational Efficiency (NCE): wilcoxon test:\n")
print(NCE_wilcox)
cat("\nNormalized Computational Efficiency (NCE): Wilcox\t test:\n")
print(NCE_t_wilcox)
cat("\nNormalized Computational Efficiency (NCE): perm test:\n")
print(NCE_perm)

# Print areas
cat("Area under the power curve for t-test:\n")
print(area_t)
cat("\nArea under the power curve for permutation test:\n")
print(area_wilcox)
cat("Area under the power curve for t-test/Wilcoxon:\n")
print(area_t_wilcox)
cat("\nArea under the power curve for permutation test:\n")
print(area_perm)

save(sample_sizes, error_t_wilcox_table, error_perm_table, error_t_table, error_wilcox_table,
                       time_t_table, time_wilcox_table, time_t_wilcox_table, time_perm_table,
                        time_diff_t, time_diff_wilcox, time_diff_perm,
                        NCE_t, NCE_wilcox, NCE_t_wilcox, NCE_perm,
                       area_t, area_wilcox_perm, area_t_wilcox, area_perm,
                      file = "TwoSampleCompareTestMethods_error.RData" )
