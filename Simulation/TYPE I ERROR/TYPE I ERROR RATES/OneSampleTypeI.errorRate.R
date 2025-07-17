#rm(list = ls())
# # set directories in cluster
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

# set directories in local computer
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("/Users/benedictkongyir/Library/Mobile Documents/com~apple~CloudDocs/PhD Thesis/Type I error")

# set up cores for parallel processing
par_set <- function(cores_reserve = 2) 
{
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if(Sys.info()["sysname"] == "Windows"){
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)     
    
  }else{
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

## Set up the simulation parameters
{
  Nsim <- 1e4;  
  alpha <- 0.05
  distributions <- c("Normal", "Uniform", "t", "Contaminated","Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal",  "Pareto")
  sample_size <- c(8, 10, 15, 20, 25, 30, 50)
}

# Progress taskbar setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size)* length(distributions) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
## Perform simulation
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"), 
                     .options.snow=opts) %:%
    foreach(dist = distributions)%dopar%
    {
      set.seed(12345) 
      pval_t <- pval_wilcox <- numeric(Nsim) 
      
      for (i in 1:Nsim) {
        x <- generate_data(n, dist) 
        pval_t[i] <- t.test(x)$p.value
        pval_wilcox[i] <- wilcox.test(x)$p.value
      }
      
      error_t_test    <- mean(pval_t < alpha)
      error_wilcox_test    <- mean(pval_wilcox < alpha)
      
      results <- list(
            error_t_test = error_t_test,
            error_wilcox_test = error_wilcox_test
        )
    }
  
  close_cluster(my_cl)        
  
  # Organize Results in tables
  errorvec <- numeric(length(sample_size) * length(distributions))
  TypeI_error_t.test <- TypeI_error_wilcox.test <- array(errorvec, 
           dim = c(length(sample_size),length(distributions)), dimnames = list(sample_size, distributions))
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      
        TypeI_error_t.test[i, j] <- (sim_out[[i]][[j]]$error_t_test)
        TypeI_error_wilcox.test[i, j] <- (sim_out[[i]][[j]]$error_wilcox_test)
      
    }
  }
})
# print results
cat("Type I error Rates for t-test")
print(TypeI_error_t.test)

cat("Type I error Rates for Wilcoxon-test")
print(TypeI_error_wilcox.test)

# save data
save(
  Nsim,
  distributions,
  sampl_size,
  TypeI_error_t.test,
  TypeI_error_wilcox.test
)

