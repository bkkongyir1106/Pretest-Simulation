#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# clear working environment
rm(list = ls())

# set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/TYPE I ERROR/TYPE I ERROR RATES")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Load the required package
if (!require(BayesFactor)) {
  install.packages("BayesFactor", dependencies = TRUE)
}
library(BayesFactor)

# set up cores for parallel processing
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores()
  cores_use <- cores - cores_reserve
  cl <- parallel::makeCluster(cores_use)
  doParallel::registerDoParallel(cl)
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
  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Exponential", "Chi-Square", "LogNormal")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
  sig_level <- 0.05
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
                     .packages = c("LaplacesDemon", "VGAM", "BayesFactor"), 
                     .options.snow=opts) %:% 
    foreach(dist = dist_sum) %dopar% {
      set.seed(12345) # Set seed for reproducibility
      pval <- numeric(N) # collect p-values
      Bayes_errors <- 0
      
      for (i in 1:N) {
        x <- generate_data(n, dist)
        pval[i] <- t.test(x)$p.value
        
        # Perform the Bayesian one-sample t-test
        bf_result <- ttestBF(x = x, mu = 0)
        
        # Extract the Bayes factor for the alternative hypothesis
        bf <- exp(bf_result@bayesFactor$bf)
        
        # Check if BF > 1 (strong evidence for the alternative hypothesis)
        if (bf > 1.5) {
          Bayes_errors <- Bayes_errors + 1  # Increment the Type I error counter
        }
      }
      
      error <- mean(pval < alpha)
      # Calculate the Type I error rate
      Bayes_error_rate <- Bayes_errors / N
      results <- list(error = error,
                      Bayes_error_rate = Bayes_error_rate)
    }
  
  close_cluster(my_cl)         
  
  # Organize Results in tables
  errorvec <- numeric(length(nvec) * length(dist_sum))
  TypeI.errorRate <- Bayes_TypeI.error_rate <- array(errorvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
  
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      TypeI.errorRate[t, j] <- (sim_out[[t]][[j]]$error)
      Bayes_TypeI.error_rate[t, j] <- (sim_out[[t]][[j]]$Bayes_error_rate)
    }
  }
})

# print results
TypeI.errorRate
Bayes_TypeI.error_rate

# save data
# save.image(paste0("OneSampleTypeI.errorRate",".RData"))

# write results to excel
# library(writexl)
# error_dataframe <- data.frame(TypeI.errorRate)
# write_xlsx(error_dataframe, path = "OneSampleTypeI_error_rates.xlsx")
