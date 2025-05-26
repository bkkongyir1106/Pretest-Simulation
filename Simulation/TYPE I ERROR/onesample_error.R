#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%-----
# clear working environment
#rm(list = ls())
# #set directories in cluster
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")

# # set directories in local computer
# setwd("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/TYPE I ERROR")
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


## Set up the simulation parameters
{
  N <- 1e4  
  alpha <- 0.05
  distribution <- c("Normal", "Uniform", "Exponential", "LogNormal")
  sample_size <- c(8, 10, 15, 20, 25, 30, 40)
}

# Progress taskbar setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size) * length(distribution) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
## Perform simulation
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow=opts, .options.RNG = 12345) %:%
            foreach(dist = distribution) %dopar%{
              pval<- numeric(N) 
              for (i in 1:N) {
                x <- generate_data(n, dist) 
                pval[i] <- t.test(x)$p.value
              }
              error    <- mean(pval < alpha)
              results <- list(error = error)
            }
          
  # Organize Results in tables
  errorvec <- numeric(length(sample_size) * length(distribution))
  OneSampleTypeI.errorRate_u <- array(errorvec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
  for (t in seq_along(sample_size)) {
    for (j in seq_along(distribution)) {
      
        OneSampleTypeI.errorRate_u[t, j] <- (sim_out[[t]][[j]]$error)

    }
  }
})
# print results
OneSampleTypeI.errorRate_u

# save data
#save(sample_size, N, OneSampleTypeI.errorRate_u, file = "OneSampleTypeI.errorRate_u.RData")

# ------------------- One Sample Two-step Type I error -----------------------
## Perform simulation
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow=opts, .options.RNG = 12345) %:%
    foreach(dist = distribution) %dopar%{
      
      pval<- numeric(N)
      for(i in 1: N) {
        x <- generate_data(n, dist) 
        if(shapiro.test(x)$p.value > alpha)
        {
          pval[i] <- t.test(x)$p.value
        }else{
          pval[i] <- wilcox.test(x)$p.value
        }
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  
  ## Organize Results in tables
  errorvec <- numeric(length(sample_size) * length(distribution) )
  OneSampleTypeI.errorRate_2step <- array(errorvec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
  for (t in seq_along(sample_size)) {
    for (j in seq_along(distribution)) {
      
        OneSampleTypeI.errorRate_2step[t, j] <- (sim_out[[t]][[j]]$error)
      
    }
  }
})

# print results
OneSampleTypeI.errorRate_2step
# save data
#save(sample_size, N, OneSampleTypeI.errorRate_2step, file = "OneSampleTypeI.errorRate_2step.RData")


# -------------- One Sample Conditional Type I error-----------------------
## Perform simulation
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow=opts, .options.RNG = 12345) %:%
    foreach(dist = distribution) %dopar%{ 
      pval<- numeric(N)  
      SamplePast = 0 
      TotalSim = 0   
      while (SamplePast < N) {
        x <- generate_data(n, dist) 
        TotalSim = TotalSim + 1
        if(shapiro.test(x)$p.value > alpha)
        {
          SamplePast = SamplePast + 1
          pval[SamplePast] <- t.test(x)$p.value
        }
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  
  # Organize Results in tables
  errorvec <- numeric(length(sample_size) * length(distribution) )
  OneSampleTypeI.errorRate_c <- array(errorvec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
  for (t in seq_along(sample_size)) {
    for (j in seq_along(distribution)) {
    
        OneSampleTypeI.errorRate_c[t, j] <- (sim_out[[t]][[j]]$error)
      
    }
  }
})
# print results
OneSampleTypeI.errorRate_c

# save data
#save(sample_size, N, OneSampleTypeI.errorRate_c, file = "OneSampleTypeI.errorRate_c.RData")

close_cluster(my_cl) 
