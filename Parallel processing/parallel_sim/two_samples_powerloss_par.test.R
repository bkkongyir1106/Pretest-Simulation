#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%----
rm(list = ls())
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/User_defined_functions.R")
source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/utility.R")
# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
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
{
  N <- 1e3; P <- 1e4; alpha <- 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Laplace", 
                "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated")
  nvec <- c(8, 10, 15, 20, 30, 50) 
  d.vec <- c(0.25, 0.5)
  #set.seed(33)
}
# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(d.vec) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
## Perform simulation
system.time({
sim_out <- foreach(i = nvec,
                   .packages = c("LaplacesDemon"),
                   .options.snow=opts) %:%
  foreach(dist = dist_sum) %:%
  foreach(d = d.vec) %dopar%
  {
    pval <- pval_perm <- numeric(N)
    for (k in 1:N) {
      x <- generate_data(i, dist) 
      y <- generate_data(i, dist) + d
      permuted_data <- c(x, y)
      pval[k] <- t.test(x, y, var.equal = F)$p.value
      
      observed_statistic <- calculate_test_statistic(x, y)
      # permutation test
      permuted_statistics <- rep(0, P)
      for (l in 1:P) 
      {
        sampled_data <- sample(permuted_data)
        permuted_data1 <- sampled_data[1:length(x)]
        permuted_data2 <- sampled_data[(length(x) + 1):(length(x) + length(y))]
        permuted_statistics[l] <- calculate_test_statistic(permuted_data1, permuted_data2)
        
      }
      pval_perm[k] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
    }
    
    power_t_test      <- mean(pval < alpha)
    power_permutation <- mean(pval_perm < alpha)  
    Results <- list(
      power_t_test= power_t_test,
      power_permutation = power_permutation
    )
  }
})
close_cluster(my_cl)        

## Output
powervec <- numeric(length(nvec) * length(dist_sum) * length(d.vec))
power_t_test_par <- power_permutation_par <- array(powervec,dim = c(length(nvec), length(dist_sum), 
                                              length(d.vec)), dimnames = list(nvec, dist_sum, d.vec))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (i in seq_along(d.vec)) {
      power_t_test_par[t, j, i] <- (sim_out[[t]][[j]][[i]]$power_t_test)
      power_permutation_par[t, j, i] <- (sim_out[[t]][[j]][[i]]$power_permutation)
    }
  }
}
power_loss_par <- power_permutation_par - power_t_test_par
print(power_t_test_par)
print(power_permutation_par)
print(power_loss_par)

save.image(paste0("TwoSamples_powerloss_par",".RData"))
