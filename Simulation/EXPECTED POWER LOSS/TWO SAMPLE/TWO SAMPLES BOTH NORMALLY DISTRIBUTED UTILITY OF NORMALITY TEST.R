#rm(list = ls())
# setwd("/home/kongyir/spring2024/power")
# source("/home/kongyir/spring2024/User_defined_functions.R")
# source("/home/kongyir/spring2024/utility.R")
setwd("~/Desktop/OSU/Research/Literature/plots")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")


# =============== Parallel process setup =========----
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
  P <- 1e3 
  alpha <- 0.05
  distribution <- c("Normal", "Uniform", "Exponential", "LogNormal")
  #distribution = "Exponential"
  sample_sizes <- c(5, 10, 15, 20, 25, 30, 40, 50) 
  effect_size <-  0.5
}

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_sizes) * length(distribution) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}
# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

system.time({
  sim_out <- foreach(n = sample_sizes,
                     .packages = c("LaplacesDemon", "VGAM"), 
                     .options.snow=opts) %:% 
            foreach(dist = distribution) %dopar%
    {
      set.seed(1234) 
      TotalSim.passed.SW.test = 0  
      TotalSim = 0 
      pval_powr = pval_perm_powr  = c()
      pval_error = pval_perm_error  = c()
      while (TotalSim.passed.SW.test < N) {
        x <- generate_data(n, dist) 
        y <- generate_data(n, dist)
        TotalSim <- TotalSim + 1 
        if(shapiro.test(x)$p.value < alpha & shapiro.test(y)$p.value < alpha){
          TotalSim.passed.SW.test <- TotalSim.passed.SW.test + 1
          # t test
          pval_powr[TotalSim.passed.SW.test] <- t.test(x, y + effect_size)$p.value
          pval_error[TotalSim.passed.SW.test] <- t.test(x, y )$p.value
          # Permutation power test
          observed_statistic <- calculate_test_statistic(x, y + effect_size)
          permuted_statistics <- rep(0, P)
          for (l in 1 : P) { 
            sampled_data <- sample(c(x, y + effect_size))
            permuted_data1 <- sampled_data[1:length(x)]
            permuted_data2 <- sampled_data[(length(x) + 1):(length(x) + length(y))]
            permuted_statistics[l] <- calculate_test_statistic(permuted_data1, permuted_data2)
          }
          pval_perm_powr[TotalSim.passed.SW.test] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
          
          # Permutation error test
          observed_statistic_error <- calculate_test_statistic(x, y)
          permuted_statistics_error <- rep(0, P)
          for (l in 1 : P) { 
            sampled_data_error <- sample(c(x, y ))
            permuted_data1_error <- sampled_data_error[1:length(x)]
            permuted_data2_error <- sampled_data_error[(length(x) + 1):(length(x) + length(y))]
            permuted_statistics_error[l] <- calculate_test_statistic(permuted_data1_error, permuted_data2_error)
          }
          pval_perm_error[TotalSim.passed.SW.test] <- round(mean(abs(permuted_statistics_error) >= abs(observed_statistic_error)), 5)
        }
      }
      power_t_test      <- mean(pval_powr < alpha)
      power_perm.test <- mean(pval_perm_powr < alpha)  
      Prob_sw <-   TotalSim.passed.SW.test/TotalSim
      
      error_t_test      <- mean(pval_error < alpha)
      error_perm.test <- mean(pval_perm_error < alpha)  
      Results <- list(
        power_t_test = power_t_test,
        power_perm.test = power_perm.test,
        error_t_test = error_t_test,
        error_perm.test = error_perm.test,
        Prob_sw = Prob_sw
      )
      return(Results)
    }
})
close_cluster(my_cl)        

# =========== Output Results ====================--

power_perm<- power_t <- p_sw <- array(NA, dim = c(length(sample_sizes), length(distribution)), dimnames = list(sample_sizes, distribution))
error_perm<- error_t <- array(NA, dim = c(length(sample_sizes), length(distribution)), dimnames = list(sample_sizes, distribution))
for (t in seq_along(sample_sizes)) {
  for (j in seq_along(distribution)) {
      power_t[t, j] <- (sim_out[[t]][[j]]$power_t_test)
      power_perm[t, j] <- (sim_out[[t]][[j]]$power_perm.test)
      error_t[t, j] <- (sim_out[[t]][[j]]$error_t_test)
      error_perm[t, j] <- (sim_out[[t]][[j]]$error_perm.test)
      p_sw[t, j] <- round((sim_out[[t]][[j]]$Prob_sw), 4)
  }
}
# print results
# power
print("Power loss")
powerloss <- power_perm - power_t
power_t
power_perm
powerloss
p_sw
ExpectedPowerloss = powerloss * p_sw
ExpectedPowerloss

# Type I error
print("Type I error Inflation")
inflation <- error_t - alpha
error_t
error_perm
inflation
Expectedinflation = inflation * p_sw
Expectedinflation

save.image(paste0("TwoSamplesExpectedPowerloss_error_04152025",".RData"))
