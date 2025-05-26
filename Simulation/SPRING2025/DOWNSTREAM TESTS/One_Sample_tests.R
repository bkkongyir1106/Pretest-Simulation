rm(list = ls())
# set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS")
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

#--------------------------------------------------------------
#  Set up the simulation parameters
# ------------------------------------------------------------
{
  N <- 1e1;  alpha <- 0.05
  #dist_sum <- c("Standard Normal", "Uniform", "t", "Contaminated","Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal",  "Pareto")
  dist_sum <- c("Standard Normal", "Uniform", "t", "Contaminated","Laplace","Exponential", "Chi-Square", "Gamma", "Weibull", "LogNormal")
  nvec <- c(8, 10, 15, 20, 25, 30) #, 50)
  sig_level <- c(0.05)
  effect_size <- 0.5
}

# ---------------------------------------------------------
# Progress taskbar setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(sig_level) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
# -------------- One Sample Unconditional t-test --------------
system.time({
  sim_out <- foreach(n = nvec,
                     .packages = c("LaplacesDemon", "VGAM"), 
                     .options.snow=opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(alpha = sig_level) %dopar%
    {
      set.seed(12345) 
      pval<- numeric(N) 
      for (i in 1:N) {
        x <- generate_data(n, dist) 
        pval[i] <- t.test(x + effect_size)$p.value
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  # Organize Results in tables
  errorvec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
  OneSampleTypeI.errorRate_u <- array(errorvec, dim = c(length(nvec), length(dist_sum), length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      for (k in seq_along(sig_level)) {
        OneSampleTypeI.errorRate_u[t, j, k] <- (sim_out[[t]][[j]][[k]]$error)
      }
    }
  }
})
# print results
OneSampleTypeI.errorRate_u

# save data
save(nvec, N, OneSampleTypeI.errorRate_u, file = "OneSampleTypeI.errorRate_u.RData")

# -------------------------------------------------------------
# ---------------One Sample Two-step t-test------ -------------
# -------------------------------------------------------------
system.time({
  sim_out <- foreach(n = nvec,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow=opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(alpha = sig_level) %dopar%
    {
      set.seed(12345) # Set seed for reproducibility
      pval<- numeric(N)
      for(i in 1: N) {
        x <- generate_data(n, dist) 
        if(shapiro.test(x)$p.value > alpha)
        {
          pval[i] <- t.test(x + effect_size)$p.value
        }else{
          pval[i] <- wilcox.test(x + effect_size)$p.value
        }
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  ## Organize Reseults in tables
  errorvec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
  OneSampleTypeI.errorRate_2step <- array(errorvec, dim = c(length(nvec), length(dist_sum), length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      for (k in seq_along(sig_level)) {
        OneSampleTypeI.errorRate_2step[t, j, k] <- (sim_out[[t]][[j]][[k]]$error)
      }
    }
  }
})

# print results
OneSampleTypeI.errorRate_2step
# save data
save(nvec, N, OneSampleTypeI.errorRate_2step, file = "OneSampleTypeI.errorRate_2step.RData")

# ----------------------------------------------------------------
# ----------------One Sample Conditional t-test------ ------------
# ----------------------------------------------------------------
system.time({
  sim_out <- foreach(n = nvec,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow=opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(alpha = sig_level) %dopar%
    {
      set.seed(12345) 
      pval<- numeric(N)  
      SamplePast = 0 # Set counter for number of samples that past SW test
      TotalSim = 0   # Set counter for total number of samples generated
      while (SamplePast < N) {
        x <- generate_data(n, dist) 
        TotalSim = TotalSim + 1
        if(shapiro.test(x)$p.value > alpha)
        {
          SamplePast = SamplePast + 1
          pval[SamplePast] <- t.test(x + effect_size)$p.value
        }
      }
      error    <- mean(pval < alpha)
      results <- list(error = error)
    }
  
  # Organize Results in tables
  errorvec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
  OneSampleTypeI.errorRate_c <- array(errorvec, dim = c(length(nvec), length(dist_sum), length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      for (k in seq_along(sig_level)) {
        OneSampleTypeI.errorRate_c[t, j, k] <- (sim_out[[t]][[j]][[k]]$error)
      }
    }
  }
})
# print results
OneSampleTypeI.errorRate_c

# save data
save(nvec, N, OneSampleTypeI.errorRate_c, file = "OneSampleTypeI.errorRate_c.RData")

close_cluster(my_cl) 

