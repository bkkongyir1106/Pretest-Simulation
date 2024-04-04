#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
rm(list = ls())
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/User_defined_functions.R")
source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/utility.R")

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
  dist_sum <- c("Standard Normal", "Uniform", "t", "Laplace","Exponential", 
                "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated")
  testvec <- c("KS", "SW", "JB", "DAP", "AD", "SF")
  nvec <- c(8, 10, 15, 20, 30, 50)
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
                   .packages = c("LaplacesDemon", "nortest", "dgof", "moments"),
                   .options.snow=opts) %:%
                foreach(dist = dist_sum) %:%
                foreach(testvec= testvec) %dopar%
  {
    pval_test_x <- pval_test_y <- numeric(N)
    for (k in 1:N) {
      x <- generate_data(n, dist)  
      y <- generate_data(n, dist)
      pval_test_x[k] <- generate_tests(x, testvec)$p.value
      pval_test_y[k] <- generate_tests(y, testvec)$p.value
    }
    powr <- mean(pval_test_x < alpha/2 | pval_test_y <alpha/2)
    results <- list(powr = powr)
  }

close_cluster(my_cl)        

## Output
powervec <- numeric(length(nvec) * length(dist_sum) * length(testvec))
power_norm.test_par <- array(powervec, dim = c(length(nvec), 
            length(dist_sum), length(testvec)), dimnames = list(nvec, dist_sum, testvec))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (i in seq_along(testvec)) {
      power_norm.test_par[t, j, i] <- (sim_out[[t]][[j]][[i]]$powr)
    }
  }
}
})
power_norm.test_par

save.image(paste0("TwoSamples_power_norm_par",".RData"))