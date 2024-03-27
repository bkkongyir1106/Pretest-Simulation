#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%----
rm(list = ls())
# Setting seed for reproducibility
set.seed(1234)  
# source("~/Desktop/OSU/Research/Permutation Test/New Simulations/cluster/User_defined_functions.R")
# source("~/Desktop/OSU/Research/Pretest-Simulation/Permutation test/R scripts/utility.R")
source("/home/kongyir/spring2024/User_defined_functions.R")
source("/home/kongyir/spring2024/utility.R")
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
  parallel::stopCluster(cl) # close the cluster
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%----
## Set up
{
  N <- 1e4; P <- 1e4; alpha <- 0.05
  dist_sum <- c("Standard Normal",  "Uniform", "t", "Laplace", "Contaminated")
  nvec <- c(8, 10, 15, 20, 30, 50) 
  sig_level <- c(0.01, 0.05, 0.075)
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
#set.seed(33)
sim_out <- foreach(i = nvec,
                   .packages = c("LaplacesDemon"),
                   .options.snow=opts,
                   .inorder=FALSE,  # This ensures that seeds are independent
                   .verbose=FALSE) %:%
  foreach(dist = dist_sum) %:%
  foreach(alpha = sig_level) %dopar%
  {
    set.seed(NULL) # Setting seed for each parallel worker
    pval <- pval_perm <- numeric(N)
    for (k in 1:N) {
      x <- generate_data(i, dist) 
      pval[k] <- t.test(x)$p.value
      
      observed_statistic <- calculate_test_statistic(x)
      
      permuted_statistics <- rep(0, P)
      for (l in 1:P) {
            myIndex <- sample(c(-1, 1), length(x), replace = TRUE)
            sample_data <- myIndex * abs(x)
            permuted_statistics[l] <- calculate_test_statistic(sample_data)
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
error_inflation <- TypeI_error_t.test - TypeI_error_perm.test
print(TypeI_error_t.test)
print(TypeI_error_perm.test)
print(error_inflation)

save.image(paste0("OneSample_error_inflation",".RData"))

