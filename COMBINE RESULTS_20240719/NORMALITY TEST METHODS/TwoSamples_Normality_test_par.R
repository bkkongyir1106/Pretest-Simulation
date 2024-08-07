#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# clear working environment
rm(list = ls())
# set directories in cluster
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

# set directories in local computer
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# set up cores for parallel processing
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
{
  N <- 1e4;  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated", "Pareto")
  testvec <- c("KS", "SW", "JB", "DAP", "AD", "SF", "CVM")
  nvec <- c(8, 10, 15, 20, 25, 30, 50)
}

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(testvec) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
## Perform simulation
system.time({
sim_out <- foreach(n = nvec,
                   .packages = c("LaplacesDemon", "VGAM", "moments", "nortest"),
                   .options.snow=opts) %:%
                foreach(dist = dist_sum) %:%
                foreach(testvec= testvec) %dopar%
  {
    pval_test_x <- pval_test_y <- numeric(N)
    for (i in 1:N) {
      x <- generate_data(n, dist)  
      y <- generate_data(n, dist)
      pval_test_x[i] <- generate_tests(x, testvec)$p.value
      pval_test_y[i] <- generate_tests(y, testvec)$p.value
    }
    powr <- mean(pval_test_x < alpha/2 | pval_test_y <alpha/2)
    results <- list(powr = powr)
  }

close_cluster(my_cl)        

## Output
powervec <- numeric(length(nvec) * length(dist_sum) * length(testvec))
power_norm.test <- array(powervec, dim = c(length(nvec), length(dist_sum), length(testvec)), dimnames = list(nvec, dist_sum, testvec))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (k in seq_along(testvec)) {
      power_norm.test[t, j, k] <- (sim_out[[t]][[j]][[k]]$powr)
    }
  }
}
})
power_norm.test

save.image(paste0("TwoSamplePowerNormalityTest",".RData"))
