# clear working environment
rm(list = ls())
# set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/TYPE I ERROR RATES/Selective Inference Type I error")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# set directories in cluster
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

# set up cores for parallel processing
par_set <- function(cores_reserve = 2) 
{
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if(Sys.info()["sysname"] == "Windows"){
    cl <- parallel::makeCluster(cores_use) # make a socket cluster
    doParallel::registerDoParallel(cl)     # for both Windows & Unix-like
    
  }else{
    
    cl <- snow::makeSOCKcluster(cores_use) # make a socket cluster
    doSNOW::registerDoSNOW(cl)           # for  
  }
  foreach::getDoParWorkers()
  return(cl)
}

close_cluster <- function(cl) {
  parallel::stopCluster(cl) # close the cluster
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%----
## Set up
{
  N <- 1e4;  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Contaminated", "Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal",  "Pareto")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
  sig_level <- c(0.05)
}

# Parallelized simulation setup
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
                     .packages = c("LaplacesDemon", "VGAM", "moments"),
                     .options.snow=opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(alpha = sig_level) %dopar%
    {
      set.seed(12345)
      pval<- numeric(N)
      for (i in 1:N) {
        x <- generate_data(n, dist) 
        y <-generate_data(n, dist)
        pval[i] <- t.test(x, y, var.equal= F)$p.value
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  close_cluster(my_cl)        
  
  ## Output
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
TypeI.errorRate

save.image(paste0("TwoSampleTypeI.errorRate",".RData"))
# Write data to excel
library(writexl)
error_dataframe <- data.frame(TypeI.errorRate)
write_xlsx(error_dataframe, path = "TwoSampleTypeI_error_rates.xlsx")
