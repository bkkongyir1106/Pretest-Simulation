#rm(list = ls())
setwd("/home/kongyir/spring2024/power")
source("/home/kongyir/spring2024/User_defined_functions.R")
source("/home/kongyir/spring2024/utility.R")
# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research.05.10.2024/Power")
# source("~/Desktop/OSU/Research/Pretest-Simulation/User defined functions/User_defined_functions.R")
# source("~/Desktop/OSU/Research/Pretest-Simulation/User defined functions/utility.R")

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
# set vectors
N <- 1e4; P <- 1e5; alpha <- 0.05
dist_sum <- c("Standard Normal",  "Uniform", "t", "Laplace")
nvec <- c(5, 10, 15, 20, 25, 30)
d.vec <- c(0.5)

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(nvec) *  length(d.vec) * length(dist_sum) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
# Function to calculate test statistic
calculate_test_statistic <- function(x) {
  return((mean(x)*sqrt(length(x))) / sqrt(var(x)))
}
# %%%%%%%%%%%%%%%%%%%%% SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%----
system.time({
  sim_out <- foreach(n = nvec,
                     .packages = c("LaplacesDemon"),
                     .options.snow=opts) %:%
    foreach(dist = dist_sum) %:%
    foreach(d = d.vec) %dopar%
    {
      set.seed(1234)
      TotalSim.passed.SW.test = 0 ; TotalSim = 0; pval = pval_perm  = c()
      while (TotalSim.passed.SW.test < N) {
        x <- generate_data(n, dist) 
        TotalSim = TotalSim + 1 
        if(shapiro.test(x)$p.value > alpha){
          TotalSim.passed.SW.test = TotalSim.passed.SW.test + 1
          pval[TotalSim.passed.SW.test] <- t.test(x + d)$p.value
          observed_statistic <- calculate_test_statistic(x + d)
          permuted_statistics <- rep(0, P)
          for (l in 1:P) {
            myIndex <- sample(c(-1, 1), length(x), replace = TRUE)
            sample_data <- myIndex * abs(x + d)
            permuted_statistics[l] <- calculate_test_statistic(sample_data)
          }
          pval_perm[TotalSim.passed.SW.test] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
        }
        
      }
      power_t_test      <- mean(pval < alpha)
      power_perm.test <- mean(pval_perm < alpha)  
      Prob.SW_n.s <-   TotalSim.passed.SW.test/TotalSim
      Results <- list(
        power_t_test = power_t_test,
        power_perm.test = power_perm.test,
        Prob.SW_n.s = TotalSim.passed.SW.test/TotalSim
      )
    }
})
close_cluster(my_cl)        

# %%%%%%%%%%%%%%%%%%%%%%%%%%%% Output Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%-----
powervec <- numeric(length(nvec) * length(dist_sum) * length(d.vec))
power_perm.test_par <- power_t.test_par <- prob.non.sig.SW.test_par <- array(powervec,dim = c(length(nvec), length(dist_sum), 
          length(d.vec)), dimnames = list(nvec, dist_sum, d.vec))
for (t in seq_along(nvec)) {
  for (j in seq_along(dist_sum)) {
    for (i in seq_along(d.vec)) {
      power_t.test_par[t, j, i] <- (sim_out[[t]][[j]][[i]]$power_t_test)
      power_perm.test_par[t, j, i] <- (sim_out[[t]][[j]][[i]]$power_perm.test)
      prob.non.sig.SW.test_par[t, j, i] <- (sim_out[[t]][[j]][[i]]$Prob.SW_n.s)
    }
  }
}
powerloss <- power_perm.test_par - power_t.test_par
power_t.test_par
power_perm.test_par
powerloss
prob.non.sig.SW.test_par

## Save results
save.image(paste0("OneSample_Expected_powerloss_par_05_10_2024",".RData"))