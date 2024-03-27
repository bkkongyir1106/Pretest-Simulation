#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%----
rm(list = ls())
# Setting seed for reproducibility
set.seed(1234)  
# source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/User_defined_functions.R")
# source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/parallel_sim/utility.R")
source("/home/kongyir/spring2024/User_defined_functions.R")
source("/home/kongyir/spring2024/utility.R")
# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

par_set <- function(cores_reserve = 2, seed = NULL) 
{
    if(!is.null(seed)){
      set.seed(seed) # setting the seed
    }
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if(Sys.info()["sysname"] == "Windows"){
    cl <- parallel::makeCluster(cores_use) # make a socket cluster
    doParallel::registerDoParallel(cl)     # for both Windows & Unix-like
    
  }else{
    
    cl <- snow::makeSOCKcluster(cores_use) # make a socket cluster
    doSNOW::registerDoSNOW(cl)           # for  
    
    #OR
    #doMC::registerDoMC(cores = cores_use)   # make a fork cluster
  }
  foreach::getDoParWorkers()
  return(cl)
}

close_cluster <- function(cl) {
  parallel::stopCluster(cl) # function to close the cluster
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%----
## Set up
{
  N <- 1e4; P <- 1e4 
  dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Laplace", 
               "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated")
  nvec <- c(8, 10, 15, 20, 30, 50) 
  sig_level <- c(0.01, 0.05, 0.075)
  #set.seed(33)
}

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2, seed = 1234)
  ntasks <- length(nvec) *  length(sig_level) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
}
## Perform simulation
system.time({
sim_out <- foreach(nvec = nvec, 
                   .packages = c("LaplacesDemon"),
                   .options.snow=opts) %:%
                  foreach(dist = dist_sum) %:%
                  foreach(alpha = sig_level) %dopar%
  {
    pval <- pval_perm <- numeric(N)
    for (k in 1:N) {
      x <- generate_data(nvec, dist) 
      y <- generate_data(nvec, dist) 
      data <- c(x, y)
      pval[k] <- t.test(x, y, var.equal = F)$p.value
      
      observed_statistic <- calculate_test_statistic(x, y)
      
      permuted_statistics <- rep(0, P)
      for (l in 1:P) 
      {
        sampled_data <- sample(data)
        permuted_data1 <- sampled_data[1:length(x)]
        permuted_data2 <- sampled_data[(length(x) + 1):(length(x) + length(y))]
        permuted_statistics[l] <- calculate_test_statistic(permuted_data1, permuted_data2)
        
      }
      pval_perm[k] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
    }
    
    error_t.test      <- mean(pval < alpha)
    error_perm.test <- mean(pval_perm < alpha)  
    Results <- list(
      error_t.test = error_t.test,
      error_perm.test = error_perm.test
    )
  }
close_cluster(my_cl)        

## Output
powervec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
TypeI_error_t.test <- TypeI_error_perm.test <- array(powervec,dim = c(length(nvec), length(dist_sum), 
                                 length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (i in seq_along(sig_level)) {
      TypeI_error_t.test[t, j, i] <- (sim_out[[t]][[j]][[i]]$error_t.test)
      TypeI_error_perm.test[t, j, i] <- (sim_out[[t]][[j]][[i]]$error_perm.test)
    }
  }
}
})
error_inflation <- TypeI_error_t.test - TypeI_error_perm.test
print(TypeI_error_t.test)
print(TypeI_error_perm.test)
print(error_inflation)

save.image(paste0("TwoSamples_error_inflation_par",".RData"))

