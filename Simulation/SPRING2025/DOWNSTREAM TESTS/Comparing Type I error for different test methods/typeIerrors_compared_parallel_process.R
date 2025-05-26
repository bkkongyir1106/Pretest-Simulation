#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%----- 

# Clear the working environment
rm(list = ls())

# Load user-defined functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Set up cores for parallel processing
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores() - cores_reserve
  cl <- if (Sys.info()["sysname"] == "Windows") {
    parallel::makeCluster(cores)
  } else {
    snow::makeSOCKcluster(cores)
  }
  # Register the parallel backend
  doParallel::registerDoParallel(cl)
  return(cl)
}

# Close the cluster after the job is done
close_cluster <- function(cl) {
  parallel::stopCluster(cl)  # close the cluster
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%---- 

# Set up simulation parameters
N <- 1e4
alpha <- 0.05
distribution <- c("Normal", "Exponential", "LogNormal")
testvec <- "SW"
sample_size <- c(8, 10, 15, 20, 25, 30, 50)

# Parallelized simulation setup
my_cl <- par_set(cores_reserve = 2)
ntasks <- length(sample_size) * length(testvec) * length(distribution)
pb <- txtProgressBar(max = ntasks, style = 3)

opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# Perform simulation with parallel processing
system.time({
  sim_out <- foreach(n = sample_size, 
                     .packages = c("LaplacesDemon", "VGAM", "moments", "nortest"), 
                     .options.snow = opts) %:%
    foreach(dist = distribution) %:%
    foreach(test = testvec) %dopar% {
      set.seed(12345)  # Reset RNG for each (n, dist, test) combination
      pval <- numeric(N)
      
      # Perform N simulations and store p-values
      for (i in seq_len(N)) {
        x <- generate_data(n, dist)
        pval[i] <- generate_tests(x, test)$p.value
      }
      
      # Calculate power (proportion of p-values less than alpha)
      powr <- mean(pval < alpha)
      list(powr = powr)
    }
  
  # Close the parallel cluster
  close_cluster(my_cl)
  
  # Store results into a more efficient structure
  power_norm.test <- array(0, 
                           dim = c(length(sample_size), length(distribution), length(testvec)), 
                           dimnames = list(sample_size, distribution, testvec))
  
  # Populate the power results
  for (t in seq_along(sample_size)) {
    for (j in seq_along(distribution)) {
      for (k in seq_along(testvec)) {
        power_norm.test[t, j, k] <- sim_out[[t]][[j]][[k]]$powr
      }
    }
  }
})

# Output the results
print(power_norm.test)

# Save the results to a file
save.image("sw.test.RData")

# ==============================================================================


# Load required libraries
library(parallel)

# Load user-defined functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/Comparing Type I error for different test methods")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

# Set seed and simulation parameters
set.seed(12345)
Nsim <- 1e3
sample_size <- c(8, 10, 15, 20, 25, 30, 50)
distribution <- c("Normal", "Exponential", "LogNormal")
alpha <- 0.05
n_boot <- 1e4
effect_size <- 0.0  # Type I error = Null is true

# Initialize result matrices
t_test_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
wilcox_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
split_mat      <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
perm_mat       <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
boot_mat       <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
t_perm_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
t_boot_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))

rownames(t_test_mat) <- sample_size
colnames(t_test_mat) <- distribution
rownames(wilcox_mat) <- rownames(split_mat) <- rownames(perm_mat) <- rownames(boot_mat) <- rownames(t_perm_mat) <- rownames(t_boot_mat) <- sample_size
colnames(wilcox_mat) <- colnames(split_mat) <- colnames(perm_mat) <- colnames(boot_mat) <- colnames(t_perm_mat) <- colnames(t_boot_mat) <- distribution

# Set up a progress bar
total_iter <- length(distribution) * length(sample_size)
current_iter <- 0
pb <- txtProgressBar(min = 0, max = total_iter, style = 3)

# Start simulation
for (d in seq_along(distribution)) {
  dist_current <- distribution[d]
  
  for (s in seq_along(sample_size)) {
    n <- sample_size[s]
    
    # Run simulations in parallel
    results <- mclapply(1:Nsim, function(i) {
      set.seed(12345 + i)  # Ensure unique seed for each iteration
      x <- generate_data(n, dist_current)
      y <- generate_data(n, dist_current)
      
      # Compute p-values for all tests
      pval.t_test <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      pval.wilcox.test <- TwoSample.test(x, y, test = "Wilcox", alpha = alpha, effect_size = effect_size)
      pval.split.test <- TwoSample.test(x, y, test = "2stage.split", alpha = alpha, effect_size = effect_size)
      pval.perm.test <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      pval.boot.test <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
      
      # Adaptive t/permutation
      if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
        pval_t_perm <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      } else {
        pval_t_perm <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      }
      
      # Adaptive t/bootstrap
      if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
        pval_t_boot <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      } else {
        pval_t_boot <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
      }
      
      # Return results as a list
      list(
        t_test = pval.t_test,
        wilcox = pval.wilcox.test,
        split = pval.split.test,
        perm = pval.perm.test,
        boot = pval.boot.test,
        t_perm = pval_t_perm,
        t_boot = pval_t_boot
      )
    }, mc.cores = detectCores() - 1)  # Adjust cores as needed
    
    # Extract p-values and compute Type I error rates
    t_test_mat[s, d]   <- mean(sapply(results, function(res) res$t_test) < alpha)
    wilcox_mat[s, d]   <- mean(sapply(results, function(res) res$wilcox) < alpha)
    split_mat[s, d]    <- mean(sapply(results, function(res) res$split) < alpha)
    perm_mat[s, d]     <- mean(sapply(results, function(res) res$perm) < alpha)
    boot_mat[s, d]     <- mean(sapply(results, function(res) res$boot) < alpha)
    t_perm_mat[s, d]   <- mean(sapply(results, function(res) res$t_perm) < alpha)
    t_boot_mat[s, d]   <- mean(sapply(results, function(res) res$t_boot) < alpha)
    
    # Update progress bar
    current_iter <- current_iter + 1
    setTxtProgressBar(pb, current_iter)
  }
}

# Clean up the progress bar and stop the cluster
close(pb)

# Print the results (unchanged)
cat("\nType I Error - t-test:\n")
print(round(t_test_mat, 3))

cat("\nType I Error - Wilcoxon test:\n")
print(round(wilcox_mat, 3))

cat("\nType I Error - 2-stage split test:\n")
print(round(split_mat, 3))

cat("\nType I Error - Permutation test:\n")
print(round(perm_mat, 3))

cat("\nType I Error - Bootstrap test:\n")
print(round(boot_mat, 3))

cat("\nType I Error - Adaptive t/perm test:\n")
print(round(t_perm_mat, 3))

cat("\nType I Error - Adaptive t/bootstrap test:\n")
print(round(t_boot_mat, 3))

# Save the results
save(t_test_mat, wilcox_mat, split_mat, perm_mat, boot_mat, t_perm_mat, t_boot_mat,
     file = "two_sample.type_I_error_20250310.RData")



# ============================= One sample Case ===============================

# Set seed and simulation parameters
set.seed(12345)
Nsim <- 1e3
sample_size <- c(8, 10, 15, 20, 25, 30, 50)
distribution <- c("Normal", "Exponential", "LogNormal")
alpha <- 0.05
n_boot <- 1e4
effect_size <- 0.0  # Type I error = Null is true

# Initialize result matrices
t_test_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
wilcox_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
split_mat      <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
#perm_mat       <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
boot_mat       <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
#t_perm_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))
t_boot_mat     <- matrix(NA, nrow = length(sample_size), ncol = length(distribution))

rownames(t_test_mat) <- sample_size
colnames(t_test_mat) <- distribution
rownames(wilcox_mat) <- rownames(split_mat) <- rownames(boot_mat) <- rownames(t_boot_mat) <- sample_size
colnames(wilcox_mat) <- colnames(split_mat) <- colnames(boot_mat) <- colnames(t_boot_mat) <- distribution

# Set up a progress bar
total_iter <- length(distribution) * length(sample_size)
current_iter <- 0
pb <- txtProgressBar(min = 0, max = total_iter, style = 3)

# Start simulation
for (d in seq_along(distribution)) {
  dist_current <- distribution[d]
  
  for (s in seq_along(sample_size)) {
    n <- sample_size[s]
    
    # Run simulations in parallel
    results <- mclapply(1:Nsim, function(i) {
      set.seed(12345 + i)  # Ensure unique seed for each iteration
      x <- generate_data(n, dist_current)
      
      # Compute p-values for all tests
      pval.t_test <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size)
      pval.wilcox.test <- OneSample.test(x, test = "Wilcox", alpha = alpha, effect_size = effect_size)
      pval.split.test <- OneSample.test(x, test = "2stage.split", alpha = alpha, effect_size = effect_size)
      #pval.perm.test <- OneSample.test(x, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      pval.boot.test <- one_sample_bootstrap_p_value(x, effect_size = effect_size, n_bootstrap = n_boot)
      
      # # Adaptive t/permutation
      # if (shapiro.test(x)$p.value > alpha) {
      #   pval_t_perm <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size)
      # } else {
      #   pval_t_perm <- OneSample.test(x, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      # }
      
      # Adaptive t/bootstrap
      if (shapiro.test(x)$p.value > alpha) {
        pval_t_boot <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size)
      } else {
        pval_t_boot <- one_sample_bootstrap_p_value(x, effect_size = effect_size, n_bootstrap = n_boot)
      }
      
      # Return results as a list
      list(
        t_test = pval.t_test,
        wilcox = pval.wilcox.test,
        split = pval.split.test,
        #perm = pval.perm.test,
        boot = pval.boot.test,
        #t_perm = pval_t_perm,
        t_boot = pval_t_boot
      )
    }, mc.cores = detectCores() - 1)  # Adjust cores as needed
    
    # Extract p-values and compute Type I error rates
    t_test_mat[s, d]   <- mean(sapply(results, function(res) res$t_test) < alpha)
    wilcox_mat[s, d]   <- mean(sapply(results, function(res) res$wilcox) < alpha)
    split_mat[s, d]    <- mean(sapply(results, function(res) res$split) < alpha)
    #perm_mat[s, d]     <- mean(sapply(results, function(res) res$perm) < alpha)
    boot_mat[s, d]     <- mean(sapply(results, function(res) res$boot) < alpha)
    #t_perm_mat[s, d]   <- mean(sapply(results, function(res) res$t_perm) < alpha)
    t_boot_mat[s, d]   <- mean(sapply(results, function(res) res$t_boot) < alpha)
    
    # Update progress bar
    current_iter <- current_iter + 1
    setTxtProgressBar(pb, current_iter)
  }
}

# Clean up the progress bar and stop the cluster
close(pb)

# Print the results (unchanged)
cat("\nType I Error - t-test:\n")
print(round(t_test_mat, 3))

cat("\nType I Error - Wilcoxon test:\n")
print(round(wilcox_mat, 3))

cat("\nType I Error - 2-stage split test:\n")
print(round(split_mat, 3))

# cat("\nType I Error - Permutation test:\n")
# print(round(perm_mat, 3))

cat("\nType I Error - Bootstrap test:\n")
print(round(boot_mat, 3))

# cat("\nType I Error - Adaptive t/perm test:\n")
# print(round(t_perm_mat, 3))

cat("\nType I Error - Adaptive t/bootstrap test:\n")
print(round(t_boot_mat, 3))

# Save the results
save(t_test_mat, wilcox_mat, split_mat, boot_mat, t_boot_mat,
     file = "one_sample.type_I_error_20250310.RData")


## plots
# Load necessary libraries for plotting
library(ggplot2)
library(gridExtra)

# Create a data frame for Type I error rates
type1_error_df <- data.frame(
  Sample_Size = rep(sample_size, times = length(distribution)),
  Distribution = rep(distribution, each = length(sample_size)),
  t_test = as.vector(t_test_mat),
  split = as.vector(split_mat),
  perm = as.vector(perm_mat),
  t_perm = as.vector(t_perm_mat)
)

# Create a data frame for power values (Shapiro-Wilk test)
power_df <- data.frame(
  Sample_Size = rep(sample_size, times = length(distribution)),
  Distribution = rep(distribution, each = length(sample_size)),
  Power = as.vector(power_norm.test[, , "SW"])
)

# Plot Type I Errors
type1_error_plot <- ggplot(type1_error_df, aes(x = Sample_Size, y = t_test, color = Distribution)) +
  geom_line() + geom_point() +
  ggtitle("Type I Error - t-test") +
  theme_minimal() + ylab("Type I Error") + xlab("Sample Size")

split_plot <- ggplot(type1_error_df, aes(x = Sample_Size, y = split, color = Distribution)) +
  geom_line() + geom_point() +
  ggtitle("Type I Error - 2-stage Split") +
  theme_minimal() + ylab("Type I Error") + xlab("Sample Size")

perm_plot <- ggplot(type1_error_df, aes(x = Sample_Size, y = perm, color = Distribution)) +
  geom_line() + geom_point() +
  ggtitle("Type I Error - Permutation Test") +
  theme_minimal() + ylab("Type I Error") + xlab("Sample Size")

t_perm_plot <- ggplot(type1_error_df, aes(x = Sample_Size, y = t_perm, color = Distribution)) +
  geom_line() + geom_point() +
  ggtitle("Type I Error - Adaptive t/perm") +
  theme_minimal() + ylab("Type I Error") + xlab("Sample Size")

# Plot Power for SW test
power_plot <- ggplot(power_df, aes(x = Sample_Size, y = Power, color = Distribution)) +
  geom_line() + geom_point() +
  ggtitle("Power - Shapiro-Wilk Test") +
  theme_minimal() + ylab("Power") + xlab("Sample Size")

# Arrange all plots into a 2x2 grid
grid.arrange(type1_error_plot, split_plot, perm_plot, t_perm_plot, power_plot,
             ncol = 2, nrow = 3)















