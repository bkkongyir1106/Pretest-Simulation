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
  N <- 1e4
  B <- 1e3
  alpha <- 0.05
  distribution <- c("Normal", "Exponential", "Chi-Square", "LogNormal")
  #distribution <- c("LogNormal")
  sample_size <- c(8, 10, 15, 20, 25, 30, 50)
  effect_size <- 0.5
}


# Function to calculate the test statistic 
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# Function to compute the area under the curve using the trapezoidal rule
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))/(max(sample_size) - min(sample_size))
}


# Progress taskbar setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size) * length(distribution) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}

## Perform simulation
system.time({
  sim_out <- foreach(n = sample_size, .packages = c("LaplacesDemon", "VGAM"), .options.snow = opts) %:%
    foreach(dist = distribution) %dopar% {
      set.seed(12345) 
      pval_t <- pval_wilcox <- pvals <- pval_perm <- numeric(N)
      for (i in 1:N) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        # Perform t-test
          pval_t[i] <- t.test(x, y + effect_size)$p.value
        
        # Perform Wilcoxon test
          pval_wilcox[i] <- wilcox.test(x, y + effect_size)$p.value
        
        N_t = 0
          if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
            pvals[i] <- t.test(x, y + effect_size)$p.value
            N_t = N_t + 1
          } else {
            pvals[i] <- wilcox.test(x, y + effect_size)$p.value
          }
      
        # Perform permutation test
          data <- c(x, y + effect_size)
          observe_stat <- TwoSample_test_statistic(x, y + effect_size)
          permuted_stat <- numeric(B)
          for (j in 1:B) {
            sample_data <- sample(data)
            sample_x <- sample_data[1:length(x)]
            sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
            permuted_stat[j] <- TwoSample_test_statistic(sample_x, sample_y) 
          }
          pval_perm[i] <- mean(abs(permuted_stat) >= abs(observe_stat))
      }
      
      # Calculate power
      powr_t <- mean(pval_t < alpha)
      powr_wilcox <- mean(pval_wilcox < alpha)
      powr_t_wilcox <- mean(pvals < alpha)
      powr_perm <- mean(pval_perm < alpha)
      
      list(
        powr_t = powr_t,
        powr_wilcox = powr_wilcox,
        powr_t_wilcox = powr_t_wilcox,
        powr_perm = powr_perm,
        p_sw = N_t/N
      )
    }
  
  close_cluster(my_cl)
})

## Output
powrvec <- numeric(length(sample_size) * length(distribution))
power_t <- power_wilcox <- power_t_wilcox <- power_perm <- p <- array(powrvec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))

for (t in seq_along(sample_size)) {
  for (j in seq_along(distribution)) {
    power_t[t, j] <- (sim_out[[t]][[j]]$powr_t)
    power_wilcox[t, j] <- (sim_out[[t]][[j]]$powr_wilcox)
    power_t_wilcox[t, j] <- (sim_out[[t]][[j]]$powr_t_wilcox)
    power_perm[t, j] <- (sim_out[[t]][[j]]$powr_perm)
  }
}

# Calculate areas under the Type I error rate curves
area_t <- apply(power_t, 2, compute_area, x = sample_size)
area_wilcox <- apply(power_wilcox, 2, compute_area, x = sample_size)
area_t_wilcox <- apply(power_t_wilcox, 2, compute_area, x = sample_size)
area_perm <- apply(power_perm, 2, compute_area, x = sample_size)

# Print results
power_t
power_wilcox
power_t_wilcox
power_perm

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



