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
  N <- 1e4;  alpha <- 0.05
  dist_sum <- c("Standard Normal", "Uniform", "t", "Contaminated", "Exponential", "Laplace", "Chi-Square", "Gamma", "Weibull", "LogNormal", "Pareto")
  #dist_sum = "Gamma"
  nvec <- c(5, 10, 15, 20, 25, 30, 50) 
  d.vec <- c( 0.5)
}
# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(d.vec) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}

# Simulation setup
system.time({
  sim_out <- foreach(n = nvec, .packages = c("LaplacesDemon", "VGAM"), 
                     .options.snow = opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(d = d.vec) %dopar% {
      set.seed(1234)
      pval_t_test = c()
      pval_mann_whitney = c()
      overall_rejections = c()
      
      for(i in 1: N) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist) + d
        
        if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
          t_test_result <- t.test(x, y)
          pval_t_test[i] <- t_test_result$p.value
          overall_rejections[i] <- t_test_result$p.value < alpha
        } else {
          pval_mann_whitney[i] <- wilcox.test(x, y)$p.value
          overall_rejections[i] <- wilcox.test(x, y)$p.value < alpha
        }
      }
      # Calculate Power of tests
      power_t_test <- mean(pval_t_test < alpha, na.rm = TRUE)
      power_mann_whitney <- mean(pval_mann_whitney < alpha, na.rm = TRUE)
      overall_power <- mean(overall_rejections, na.rm = TRUE)
      
      Results <- list(
        power_t_test = power_t_test,
        power_mann_whitney = power_mann_whitney,
        overall_power = overall_power
      )
      return(Results)
    }
})
close_cluster(my_cl)

# Output Results
powervec <- numeric(length(nvec) * length(dist_sum) * length(d.vec))
power_t_test_par <- power_mann_whitney_par <- overall_power_par <- array(powervec, 
    dim = c(length(nvec), length(dist_sum), length(d.vec)), dimnames = list(nvec, dist_sum, d.vec))

for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (k in seq_along(d.vec)) {
      power_t_test_par[t, j, k] <- round((sim_out[[t]][[j]][[k]]$power_t_test), 4)
      power_mann_whitney_par[t, j, k] <- round((sim_out[[t]][[j]][[k]]$power_mann_whitney), 4)
      overall_power_par[t, j, k] <- round((sim_out[[t]][[j]][[k]]$overall_power), 4)
    }
  }
}

cat("Power of t-test:\n")
print(power_t_test_par)
cat("\nPower of Mann-Whitney U test:\n")
print(power_mann_whitney_par)
cat("\nOverall Power of the two-stage procedure:\n")
print(overall_power_par)

# Save the results
save(nvec, power_t_test_par, power_mann_whitney_par, overall_power_par, file = "TwoSamplesOverallPower20240607.RData")