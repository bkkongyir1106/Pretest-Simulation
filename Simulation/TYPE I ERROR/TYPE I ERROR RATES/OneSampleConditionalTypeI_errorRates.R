# clear working environment
# Un-comment this before submitting code to cluster
# setwd("/home/kongyir/spring2024/error")
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

#set directories in local computer
rm(list = ls())
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/TYPE I ERROR")
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
  N <- 1e5;  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t"  , "Contaminated", "Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal", "Pareto")
  nvec <- c(8, 10, 15, 20, 25, 30)
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
                     .packages = c("LaplacesDemon", "VGAM"), # Put Difficult Packages here
                     .options.snow=opts) %:%
    foreach(dist = dist_sum)%dopar%
    {
      set.seed(12345) # Set seed for resproducibility
      pval<- numeric(N)  # collect p values
      SamplePast = 0 # Set counter for number of samples that past SW test
      TotalSim = 0   # Set counter for total number of samples generated
      while (SamplePast < N) {
        x <- generate_data(n, dist) 
        TotalSim = TotalSim + 1
        if(shapiro.test(x)$p.value > alpha)
        {
          SamplePast = SamplePast + 1
          pval[SamplePast] <- t.test(x)$p.value
        }
      }
      prob_sw <- SamplePast/TotalSim
      error1 <- mean(pval < alpha)
      error    <- mean(pval < alpha * prob_sw)
      results <- list(error = error,
                      error1 = error1,
                      prob_sw = prob_sw
                      )
    }
  
  close_cluster(my_cl)        
  
  # Organize Results in tables
  errorvec <- numeric(length(nvec) * length(dist_sum))
  TypeI.errorRate <- TypeI.errorRate1 <- ProbFail_to_Reject0h <- array(errorvec, dim = c(length(nvec), length(dist_sum)), dimnames = list(nvec, dist_sum))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
        TypeI.errorRate[t, j] <- (sim_out[[t]][[j]]$error)
        TypeI.errorRate1[t, j] <- (sim_out[[t]][[j]]$error1)
        ProbFail_to_Reject0h[t, j] <- (sim_out[[t]][[j]]$prob_sw)
    }
  }
})
# print results
TypeI.errorRate
TypeI.errorRate1
ProbFail_to_Reject0h
# save Data
save.image(paste0("OneSampleCondAjustedTypeI.errorRate",".RData"))

# Write data to excel
library(writexl)
error_dataframe <- data.frame(TypeI.errorRate)
write_xlsx(error_dataframe, path = "OneConditionalTypeI_error_rates.xlsx")

