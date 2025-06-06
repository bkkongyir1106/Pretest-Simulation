#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# clear working environment
rm(list = ls())
# # set directories in cluster
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

# set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/TYPE I ERROR RATES/Selective Inference Type I error")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# set up cores for parallel processing
par_set <- function(cores_reserve = 2) 
{
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if(Sys.info()["sysname"] == "Windows"){
    # make a socket for both Windows & Unix-like
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)     
    
  }else{
    # make a socket cluster for Unix-like
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
  N <- 1e4;  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Contaminated","Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal",  "Pareto")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
  sig_level <- c(0.05)
}

# Progress taskbar setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(sig_level) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
## Perform simulation
system.time({
  sim_out <- foreach(n = nvec,
                     .packages = c("LaplacesDemon", "VGAM"), # Put Difficult Packages here
                     .options.snow=opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(alpha = sig_level) %dopar%
    {
      set.seed(12345) # Set seed for resproducibility
      pval<- numeric(N) # collect p values
      for (i in 1:N) {
        x <- generate_data(n, dist) 
        pval[i] <- t.test(x)$p.value
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  close_cluster(my_cl)        
  
  # Organize Results in tables
  errorvec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
  TypeI.errorRate <- array(errorvec, dim = c(length(nvec), length(dist_sum), length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      for (k in seq_along(sig_level)) {
        TypeI.errorRate[t, j, k] <- (sim_out[[t]][[j]][[k]]$error)
      }
    }
  }
})
# print results
TypeI.errorRate

# save data
save.image(paste0("OneSampleTypeI.errorRate",".RData"))

# write results to excel
library(writexl)
error_dataframe <- data.frame(TypeI.errorRate)
write_xlsx(error_dataframe, path = "OneSampleTypeI_error_rates.xlsx")
