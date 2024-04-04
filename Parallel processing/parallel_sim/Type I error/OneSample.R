#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%----
rm(list = ls())
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Parallel processing/parallel_sim/Type I error")
source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/User_defined_functions.R")
source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/utility.R")
# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x) {
  return((mean(x)*sqrt(length(x))) / sqrt(var(x)))
}
par_set <- function(cores_reserve = 2) 
{
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if(Sys.info()["sysname"] == "Windows"){
    cl <- parallel::makeCluster(cores_use) # make a socket cluster
    doParallel::registerDoParallel(cl)     # for both Windows & Unix-like
    
    # OR
    #cl <- snow::makeSOCKcluster(cores_use) # make a socket cluster
    #doSNOW::registerDoSNOW(cl)           # for  Windows only    
    
  }else{
    #cl <- parallel::makeCluster(cores_use) # make a socket cluster
    #doParallel::registerDoParallel(cl)     # for Windows & Unix-like
    
    cl <- snow::makeSOCKcluster(cores_use) # make a socket cluster
    doSNOW::registerDoSNOW(cl)           # for  
    
    #OR
    #doMC::registerDoMC(cores = cores_use)   # make a fork cluster
  }
  foreach::getDoParWorkers()
  return(cl)
}

close_cluster <- function(cl) {
  parallel::stopCluster(cl) # close the cluster
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%----
## Set up
## Set up
{
  N <- 1e4; alpha <- 0.05; treshold = 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Laplace", 
                "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated")
  nvec <- c(8, 10, 15, 20, 30, 50) 
  sig_level <- c(0.01, 0.05, 0.075)
}
# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(dist_sum) * length(sig_level) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
## Perform simulation
sim_out <- foreach(nvec = nvec,
                   .packages = c("LaplacesDemon"),
                   .options.snow=opts) %:%
  foreach(dist = dist_sum) %:%
  foreach(alpha = sig_level) %dopar%
  {
    pval <- numeric(N)
    for (k in 1:N) {
      x <- generate_data(nvec, dist) 
      pval[k] <- t.test(x)$p.value
    }
    
    error_t.test      <- mean(pval < alpha)
    Results <- list(
      error_t.test = error_t.test
    )
  }
    
# Stop the cluster when done
stopCluster(my_cl)
## Output
powervec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
TypeI_error_t.test <- array(powervec,dim = c(length(nvec), length(dist_sum), 
               length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (i in seq_along(sig_level)) {
      TypeI_error_t.test[t, j, i] <- (sim_out[[t]][[j]][[i]]$error_t.test)
    }
  }
}
error_inflation <- TypeI_error_t.test - treshold
print(TypeI_error_t.test)
print(error_inflation)

save.image(paste0("OneSample_error_inflation",".RData"))
