#rm(list = ls())
# setwd("/home/kongyir/spring2024/power")
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research.05.10.2024/Power")
source("~/Desktop/OSU/Research/Pretest-Simulation/User defined functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/User defined functions/utility.R")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parallel process setup %%%%%%%%%%%%%%%%%%%----
{
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
}

{## Set up
  N <- 1e5
  dist_sum <- c("Standard Normal", "Uniform", "t", "Contaminated", "Exponential", "Laplace", "Chi-Square", "Gamma", "Weibull", "LogNormal", "Pareto")
  #dist_sum = "Gamma"
  nvec <- c(5, 10, 15, 20, 25, 30, 50) 
  sig_level <- c( 0.05)
}
# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(sig_level) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}

# Simulation setup
system.time({
  sim_out <- foreach(n = nvec, .packages = c("LaplacesDemon", "VGAM"), 
                     .options.snow = opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(alpha = sig_level) %dopar% {
      set.seed(1234)
      pval_t_test = c()
      pval_mann_whitney = c()
      overall_rejections = c()
      
      for(i in 1: N) {
        x <- generate_data(n, dist)
        if (shapiro.test(x)$p.value > alpha) {
          t_test_result <- t.test(x)
          pval_t_test[i] <- t_test_result$p.value
          overall_rejections[i] <- t_test_result$p.value
        } else {
          mann_whitney_result <- wilcox.test(x)
          pval_mann_whitney[i] <- mann_whitney_result$p.value
          overall_rejections[i] <- mann_whitney_result$p.value
        }
      }
      # Calculate Power of tests
      error_t_test <- mean(pval_t_test < alpha, na.rm = TRUE)
      error_mann_whitney <- mean(pval_mann_whitney < alpha, na.rm = TRUE)
      overall_error <- mean(overall_rejections < alpha, na.rm = TRUE)
      
      Results <- list(
        error_t_test = error_t_test,
        error_mann_whitney = error_mann_whitney,
        overall_error = overall_error
      )
      return(Results)
    }
})
close_cluster(my_cl)

# Output Results
errorvec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
TypeI.error.t.test <- TypeI.error.U.test <- overall_TypeI.error <- array(errorvec, 
 dim = c(length(nvec), length(dist_sum), length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (k in seq_along(sig_level)) {
      TypeI.error.t.test[t, j, k] <- round((sim_out[[t]][[j]][[k]]$error_t_test), 4)
      TypeI.error.U.test[t, j, k] <- round((sim_out[[t]][[j]][[k]]$error_mann_whitney), 4)
      overall_TypeI.error[t, j, k] <- round((sim_out[[t]][[j]][[k]]$overall_error), 4)
    }
  }
}

cat("TypeI error rate of t-test:\n")
print(TypeI.error.t.test)
cat("\nTypeI error rate of Mann-Whitney U test:\n")
print(TypeI.error.U.test)
cat("\nOverall TypeI error rate of the two-stage procedure:\n")
print(overall_TypeI.error)

# Save the results
save(nvec, TypeI.error.t.test, TypeI.error.U.test, overall_TypeI.error, file = "OneSampleOverallTypeI.error.rate20240607.RData")


