# clear working environment
rm(list = ls())
# set directories in local computer
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/NORMALITY TEST METHODS")
# set up cores for parallel processing
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores()
  cores_use <- cores - cores_reserve
  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  foreach::getDoParWorkers()
  return(cl)
}

close_cluster <- function(cl) {
  parallel::stopCluster(cl) # close the cluster
}

#-------------- PERFORM SIMULATIONS ----
## Set up
{
  N <- 1e4;  alpha <- 0.05
  dist_sum <- c("Uniform", "t", "Laplace", "Exponential", "Contaminated", "LogNormal")
  testvec <- c("KS", "SW", "JB", "DAP", "AD", "SF", "CVM")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
}

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) * length(testvec) * length(dist_sum)
  pb <- txtProgressBar(max = ntasks, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
}

## Perform simulation
system.time({
  sim_out <- foreach(n = nvec,
                     .packages = c("LaplacesDemon", "VGAM", "moments", "nortest"),
                     .options.snow = opts) %:% 
    foreach(dist = dist_sum) %:% 
    foreach(test = testvec) %dopar% {   
      set.seed(12345)  
      pval <- numeric(N)
      for (i in 1:N) {
        x <- generate_data(n, dist)  
        pval[i] <- generate_tests(x, test)$p.value
      }
      powr <- mean(pval < alpha)
      list(powr = powr)
    }
  
  close_cluster(my_cl)        
  
  ## Output
  powervec <- numeric(length(nvec) * length(dist_sum) * length(testvec))
  power_norm.test <- array(powervec, 
                           dim = c(length(nvec), length(dist_sum), length(testvec)), 
                           dimnames = list(nvec, dist_sum, testvec))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      for (k in seq_along(testvec)) {
        power_norm.test[t, j, k] <- sim_out[[t]][[j]][[k]]$powr
      }
    }
  }
})
power_norm.test

save.image(paste0("OneSamplePowerNormality", ".RData"))

