#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# # Un-comment this before submitting code to cluster
# setwd("/home/kongyir")
# source("/home/kongyir/User_defined_functions.R")
# source("/home/kongyir/utility.R")

# # clear working environment
rm(list = ls())
# Set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE/POWER")
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
  #dist_sum <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal")
  dist_sum <- c("LogNormal")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
  d <- 0.5
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
  sim_out <- foreach(n = nvec, .packages = c("LaplacesDemon", "VGAM"), .options.snow = opts) %:%
    foreach(dist = dist_sum) %dopar% {
      set.seed(12345) # Set seed for reproducibility
      pval_t <- pval_wilcox <- pvals <- pval_perm <- numeric(N)
      time_t <- time_wilcox <- time_t_wilcox <- time_perm <- numeric(N)
      for (i in 1:N) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        # Perform t-test
        time_t[i] <- system.time({
          pval_t[i] <- t.test(x, y + d)$p.value
        })["elapsed"]
        
        # Perform Wilcoxon test
        time_wilcox[i] <- system.time({
          pval_wilcox[i] <- wilcox.test(x, y + d)$p.value
        })["elapsed"]
        
        # Perform t-test/Wilcoxon test based on Shapiro-Wilk normality test
        N_t = 0
        time_t_wilcox[i] <- system.time({
          if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
            pvals[i] <- t.test(x, y + d)$p.value
            N_t = N_t + 1
          } else {
            pvals[i] <- wilcox.test(x, y + d)$p.value
          }
        })["elapsed"]
        
        # Perform permutation test
        time_perm[i] <- system.time({
          data <- c(x, y + d)
          observe_stat <- calculate_test_statistic(x, y + d)
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
      
      # Calculate Type I error rates
      powr_t <- mean(pval_t < alpha)
      powr_wilcox <- mean(pval_wilcox < alpha)
      powr_t_wilcox <- mean(pvals < alpha)
      powr_perm <- mean(pval_perm < alpha)
      
      # Calculate average computation times
      avg_time_t <- mean(time_t)
      avg_time_wilcox <- mean(time_wilcox)
      avg_time_t_wilcox <- mean(time_t_wilcox)
      avg_time_perm <- mean(time_perm)
      
      list(
        powr_t = powr_t,
        powr_wilcox = powr_wilcox,
        powr_t_wilcox = powr_t_wilcox,
        powr_perm = powr_perm,
        time_t = avg_time_t,
        time_wilcox = avg_time_wilcox,
        time_t_wilcox = avg_time_t_wilcox,
        time_perm = avg_time_perm,
        p_sw = N_t/N
      )
    }
  
  close_cluster(my_cl)
})

## Output
powrvec <- numeric(length(nvec) * length(dist_sum))
power_t <- power_wilcox <- power_t_wilcox <- power_perm <- p <- array(powrvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
avg_time_t <- avg_time_wilcox <- avg_time_t_wilcox <- avg_time_perm <- array(powrvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))

for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    # Probability of Type I error rates
    power_t[t, j] <- (sim_out[[t]][[j]]$powr_t)
    power_wilcox[t, j] <- (sim_out[[t]][[j]]$powr_wilcox)
    power_t_wilcox[t, j] <- (sim_out[[t]][[j]]$powr_t_wilcox)
    power_perm[t, j] <- (sim_out[[t]][[j]]$powr_perm)
    
    # Computation Time
    avg_time_t[t, j] <- (sim_out[[t]][[j]]$time_t)
    avg_time_wilcox[t, j] <- (sim_out[[t]][[j]]$time_wilcox)
    avg_time_t_wilcox[t, j] <- (sim_out[[t]][[j]]$time_t_wilcox)
    avg_time_perm[t, j] <- (sim_out[[t]][[j]]$time_perm)
    p[t, j] <- sim_out[[t]][[j]]$p_sw
  }
}

# Calculate areas under the Type I error rate curves
area_t <- apply(power_t, 2, compute_area, x = nvec)
area_wilcox <- apply(power_wilcox, 2, compute_area, x = nvec)
area_t_wilcox <- apply(power_t_wilcox, 2, compute_area, x = nvec)
area_perm <- apply(power_perm, 2, compute_area, x = nvec)

# Print results
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
area_t
area_wilcox
area_t_wilcox
area_perm

# Save Data
#save.image(paste0("TwoSamplepowerAUC",".RData"))

# Write data to excel
# library(writexl)
# power_dataframe <- data.frame(power_t, power_wilcox, power_t_wilcox, power_perm)
# write_xlsx(power_dataframe, path = "TwoSamplepower_AUC.xlsx")



