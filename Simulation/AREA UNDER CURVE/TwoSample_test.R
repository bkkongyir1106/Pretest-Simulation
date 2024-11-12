#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# # Un-comment this before submitting code to cluster
# setwd("/home/kongyir")
# source("/home/kongyir/User_defined_functions.R")
# source("/home/kongyir/utility.R")

# clear working environment
rm(list = ls())
## Set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE")
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
  N <- 1e4
  B <- 1e4
  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
  d <- 0.75
}

# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))/(max(nvec) - min(nvec))
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
  sim_out <- foreach(n = nvec, 
                     .packages = c("LaplacesDemon", "VGAM"), 
                     .options.snow = opts) %:%
    foreach(dist = dist_sum) %dopar% {
      
      set.seed(12345)  # Set seed for reproducibility
      
      # Initialize vectors to store results
      pval_t_power <- pval_wilcox_power <- pvals_power <- pval_perm_power <- numeric(N)
      pval_t_error <- pval_wilcox_error <- pvals_error <- pval_perm_error <- numeric(N)
      pval_sw <- time_t <- time_wilcox <- time_t_wilcox <- time_perm <- numeric(N)
      
      for (i in 1:N) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        
        # Perform SW test
        pval_sw[i] <- shapiro.test(x)$p.value
        
        # Perform t-test
        time_t[i] <- system.time({
          pval_t_error[i] <- t.test(x, y)$p.value
          pval_t_power[i] <- t.test(x, y + d)$p.value
        })["elapsed"]
        
        # Perform Wilcoxon test
        time_wilcox[i] <- system.time({
          pval_wilcox_error[i] <- wilcox.test(x, y, mu = 0, alternative = "two.sided", paired = F)$p.value
          pval_wilcox_power[i] <- wilcox.test(x, y + d, mu = 0, alternative = "two.sided", paired = F)$p.value
        })["elapsed"]
        
        # Perform t-test/Wilcoxon test based on Shapiro-Wilk normality test
        time_t_wilcox[i] <- system.time({
          if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
            pvals_error[i] <- t.test(x, y)$p.value
            pvals_power[i] <- t.test(x,  y + d)$p.value
          } else {
            pvals_error[i] <- wilcox.test(x, y ,mu = 0, alternative = "two.sided", paired = F)$p.value
            pvals_power[i] <- wilcox.test(x, y + d, mu = 0, alternative = "two.sided", paired = F)$p.value
          }
        })["elapsed"]
        
        # Perform permutation test
        time_perm[i] <- system.time({
          observe_stat_error <- calculate_test_statistic(x, y)
          observe_stat_power <- calculate_test_statistic(x, y + d)
          permuted_stat_error <- permuted_stat_power <- numeric(B)
          data_error <- c(x, y)
          data_power <- c(x, y + d)
          
          for (j in 1:B) {
            sample_data_error <- sample(data_error)
            sample_data_power <- sample(data_power)
            # Calculate test statistics for permutation samples
            permuted_stat_error[j] <- calculate_test_statistic(sample_data_error[1:length(x)], 
                                                               sample_data_error[(length(x) + 1):length(data_error)])
            permuted_stat_power[j] <- calculate_test_statistic(sample_data_power[1:length(x)], 
                                                               sample_data_power[(length(x) + 1):length(data_power)])
          }
          
          # Calculate p-values
          pval_perm_error[i] <- mean(abs(permuted_stat_error) >= abs(observe_stat_error))
          pval_perm_power[i] <- mean(abs(permuted_stat_power) >= abs(observe_stat_power))
        })["elapsed"]
      }
      
      # Compute metrics and average times
      list(
        powr_sw = mean(pval_sw < alpha),
        t_error = mean(pval_t_error < alpha),
        wilcox_error = mean(pval_wilcox_error < alpha),
        t_wilcox_error = mean(pvals_error < alpha),
        perm_error = mean(pval_perm_error < alpha),
        powr_t = mean(pval_t_power < alpha),
        powr_wilcox = mean(pval_wilcox_power < alpha),
        powr_t_wilcox = mean(pvals_power < alpha),
        powr_perm = mean(pval_perm_power < alpha),
        time_t = mean(time_t),
        time_wilcox = mean(time_wilcox),
        time_t_wilcox = mean(time_t_wilcox),
        time_perm = mean(time_perm)
      )
    }
})
close_cluster(my_cl)

## Output
powrvec <- numeric(length(nvec) * length(dist_sum))
power_sw <- array(powrvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
error_t <- error_wilcox <- error_t_wilcox <- error_perm <- array(powrvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
power_t <- power_wilcox <- power_t_wilcox <- power_perm <- array(powrvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
avg_time_t <- avg_time_wilcox <- avg_time_t_wilcox <- avg_time_perm <- array(powrvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))

for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    #power of sw test
    power_sw[t, j] <- (sim_out[[t]][[j]]$powr_sw)
    # Probability of Type I error rates
    error_t[t, j] <- (sim_out[[t]][[j]]$t_error)
    error_wilcox[t, j] <- (sim_out[[t]][[j]]$wilcox_error)
    error_t_wilcox[t, j] <- (sim_out[[t]][[j]]$t_wilcox_error)
    error_perm[t, j] <- (sim_out[[t]][[j]]$perm_error)
    
    # Power
    power_t[t, j] <- (sim_out[[t]][[j]]$powr_t)
    power_wilcox[t, j] <- (sim_out[[t]][[j]]$powr_wilcox)
    power_t_wilcox[t, j] <- (sim_out[[t]][[j]]$powr_t_wilcox)
    power_perm[t, j] <- (sim_out[[t]][[j]]$powr_perm)
    
    # Computation Time
    avg_time_t[t, j] <- (sim_out[[t]][[j]]$time_t)
    avg_time_wilcox[t, j] <- (sim_out[[t]][[j]]$time_wilcox)
    avg_time_t_wilcox[t, j] <- (sim_out[[t]][[j]]$time_t_wilcox)
    avg_time_perm[t, j] <- (sim_out[[t]][[j]]$time_perm)
  }
}

# Calculate areas under the Type I error rate curves
auc_error_t <- apply(error_t, 2, compute_area, x = nvec)
auc_error_wilcox <- apply(error_wilcox, 2, compute_area, x = nvec)
auc_error_t_wilcox <- apply(error_t_wilcox, 2, compute_area, x = nvec)
auc_error_perm <- apply(error_perm, 2, compute_area, x = nvec)

# Calculate areas under the Type I error rate curves
auc_power_t <- apply(power_t, 2, compute_area, x = nvec)
auc_power_wilcox <- apply(power_wilcox, 2, compute_area, x = nvec)
auc_power_t_wilcox <- apply(power_t_wilcox, 2, compute_area, x = nvec)
auc_power_perm <- apply(power_perm, 2, compute_area, x = nvec)

#power of sw test
power_sw

# Print Probability of Type I error rates
error_t
error_wilcox
error_t_wilcox
error_perm

# Print Power results
power_t
power_wilcox
power_t_wilcox
power_perm

# Computation time
avg_time_t
avg_time_wilcox
avg_time_t_wilcox
avg_time_perm

# Area under Type I error curve
auc_error_t
auc_error_wilcox
auc_error_t_wilcox
auc_error_perm

# Area under power curve
auc_power_t
auc_power_wilcox
auc_power_t_wilcox
auc_power_perm

# Save Data
save.image(paste0("TwoSamplepowerAUC",".RData"))

# # Write data to excel
# library(writexl)
# error_dataframe <- data.frame(power_t, power_wilcox, power_t_wilcox, power_perm)
# write_xlsx(error_dataframe, path = "OneSampleTypeI_error_rates.xlsx")



