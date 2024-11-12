#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# # Un-comment this before submitting code to cluster
# setwd("/home/kongyir")
# source("/home/kongyir/User_defined_functions.R")
# source("/home/kongyir/utility.R")

# # clear working environment
rm(list = ls())
# Set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE/TYPE I ERROR")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Set up cores for parallel processing
par_set <- function(cores_reserve = 2) {
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if (Sys.info()["sysname"] == "Windows") {
    # Make a socket for both Windows & Unix-like
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)
  } else {
    # Make a socket cluster for Unix-like
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  foreach::getDoParWorkers()
  return(cl)
}

# Function to close the cluster
close_cluster <- function(cl) {
  parallel::stopCluster(cl)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%----
## Set up the simulation parameters
{
  N <- 1e3
  B <- 1e3
  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
  sig_level <- c(0.05)
}

# Function to calculate test statistic for one-sample
calculate_test_statistic <- function(x) {
  return((mean(x) * sqrt(length(x))) / sd(x))
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)) / (max(nvec) - min(nvec))
}

# Function to perform the bootstrap test
perform_bootstrap_test <- function(x, B) {
  n <- length(x)
 # Calculate the observed test statistic
  observe_stat <- calculate_test_statistic(x)
  
  # Initialize vector for bootstrap test statistics
  bootstrap_stat <- numeric(B)
  
  # Generate the bootstrap replicates
  for (j in 1:B) {
    sample_data <- sample(x, size = n, replace = TRUE)
    bootstrap_stat[j] <- calculate_test_statistic(sample_data)
  }
  
  # Calculate the p-value
  pval_bootstrap <- mean(abs(bootstrap_stat) >= abs(observe_stat))
  
  return(pval_bootstrap)
}

# Progress taskbar setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}

## Perform simulation
system.time({
  sim_out <- foreach(n = nvec, .packages = c("LaplacesDemon", "VGAM"), .options.snow = opts) %:%
    foreach(dist = dist_sum) %dopar% {
      set.seed(12345) # Set seed for reproducibility
      pval_t <- pval_bootstrap <- pvals <- numeric(N)
      time_t <- time_bootstrap <- time_t_bootstrap <- numeric(N)
      
      for (i in 1:N) {
        x <- generate_data(n, dist)
        
        # Perform t-test
        time_t[i] <- system.time({
          pval_t[i] <- t.test(x)$p.value
        })["elapsed"]
        
        # Perform bootstrap test
        time_bootstrap[i] <- system.time({
          pval_bootstrap[i] <- perform_bootstrap_test(x, B)
        })["elapsed"]
        
        # Perform t-test/bootstrap test based on Shapiro-Wilk normality test
        time_t_bootstrap[i] <- system.time({
          if (shapiro.test(x)$p.value > alpha) {
            pvals[i] <- t.test(x)$p.value
          } else {
            pvals[i] <- perform_bootstrap_test(x, B)
          }
        })["elapsed"]
      }
      
      # Calculate Type I error rates
      error_t <- mean(pval_t < alpha)
      error_bootstrap <- mean(pval_bootstrap < alpha)
      error_t_bootstrap <- mean(pvals < alpha)
      
      # Calculate average computation times
      avg_time_t <- mean(time_t)
      avg_time_bootstrap <- mean(time_bootstrap)
      avg_time_t_bootstrap <- mean(time_t_bootstrap)
      
      list(
        error_t = error_t,
        error_bootstrap = error_bootstrap,
        error_t_bootstrap = error_t_bootstrap,
        time_t = avg_time_t,
        time_bootstrap = avg_time_bootstrap,
        time_t_bootstrap = avg_time_t_bootstrap
      )
    }
  
  close_cluster(my_cl)
})

## Output
errorvec <- numeric(length(nvec) * length(dist_sum))
TypeI.errorRate_t <- TypeI.errorRate_bootstrap <- TypeI.errorRate_t_bootstrap <- array(errorvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
avg_time_t <- avg_time_bootstrap <- avg_time_t_bootstrap <- array(errorvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))

for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    # Probability of Type I error rates
    TypeI.errorRate_t[t, j] <- (sim_out[[t]][[j]]$error_t)
    TypeI.errorRate_bootstrap[t, j] <- (sim_out[[t]][[j]]$error_bootstrap)
    TypeI.errorRate_t_bootstrap[t, j] <- (sim_out[[t]][[j]]$error_t_bootstrap)
    
    # Computation Time
    avg_time_t[t, j] <- (sim_out[[t]][[j]]$time_t)
    avg_time_bootstrap[t, j] <- (sim_out[[t]][[j]]$time_bootstrap)
    avg_time_t_bootstrap[t, j] <- (sim_out[[t]][[j]]$time_t_bootstrap)
  }
}

# Calculate areas under the Type I error rate curves
area_t <- apply(TypeI.errorRate_t, 2, compute_area, x = nvec)
area_bootstrap <- apply(TypeI.errorRate_bootstrap, 2, compute_area, x = nvec)
area_t_bootstrap <- apply(TypeI.errorRate_t_bootstrap, 2, compute_area, x = nvec)

# Print results
TypeI.errorRate_t
TypeI.errorRate_bootstrap
TypeI.errorRate_t_bootstrap

# Computation time
avg_time_t
avg_time_bootstrap
avg_time_t_bootstrap

# Area under Type I error curve
area_t
area_bootstrap
area_t_bootstrap

# Save Data
#save.image(paste0("OneSampleTypeI.errorRateAUC",".RData"))

# # Write data to excel
# library(writexl)
# error_dataframe <- data.frame(TypeI.errorRate_t, TypeI.errorRate_bootstrap, TypeI.errorRate_t_bootstrap)
# write_xlsx(error_dataframe, path = "OneSampleTypeI_error_rates.xlsx")
