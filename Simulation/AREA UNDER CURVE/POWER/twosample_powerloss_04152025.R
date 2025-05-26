# ========================
# Clear environment
rm(list = ls())

# Set working directory and source custom functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE/POWER")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# ========================
# Set up parallel processing
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores()
  cores_use <- max(1, cores - cores_reserve)
  
  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)
  } else {
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

# ========================
# Parameters
N <- 1e4           
alpha <- 0.05      
distribution <- c("Normal", "Uniform", "Exponential", "LogNormal")
sample_sizes <- c(8, 10, 15, 20, 25, 30, 40, 50)
effect_size <- 0.5
B <- 1e3          

# ========================
# Setup parallel 
my_cl <- par_set(cores_reserve = 2)
parallel::clusterExport(my_cl, varlist = c("TwoSample_test_statistic", "generate_data", "alpha", "effect_size", "B"))
params <- expand.grid(n = sample_sizes, dist = distribution)
ntasks <- nrow(params)

pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# ========================
# Run simulation
system.time({
  sim_out <- foreach(i = 1:nrow(params), .packages = c("LaplacesDemon", "VGAM"), .options.snow = opts) %dopar% {
    n <- params$n[i]
    dist <- params$dist[i]
    
    set.seed(12345 + i)  
    
    Nsim <- 0
    rejectH0 <- 0
    pval <- numeric(N)
    pval_perm <- numeric(N)
    
    while (rejectH0 < N) {
      x <- generate_data(n, dist)
      y <- generate_data(n, dist)
      Nsim <- Nsim + 1
      
      if (shapiro.test(x)$p.value < alpha && shapiro.test(y)$p.value < alpha) {
        rejectH0 <- rejectH0 + 1
        
        # t-test p-value
        pval[rejectH0] <- t.test(x, y + effect_size)$p.value
        
        # Observed statistic 
        observe_stat <- TwoSample_test_statistic(x, y + effect_size)
        
        # Permutation: 
        combine_data <- c(x, y + effect_size)  
        permuted_stat <- numeric(B)
        
        for (j in 1:B) {
          sample_data <- sample(combine_data)
          sample_x <- sample_data[1:length(x)]
          sample_y <- sample_data[(length(x) + 1):(2 * length(x))]
          permuted_stat[j] <- TwoSample_test_statistic(sample_x, sample_y)
        }
        
        
        pval_perm[rejectH0] <- mean(abs(permuted_stat) >= abs(observe_stat))
      }
    }
    
    list(
      powr_t = mean(pval[1:rejectH0] < alpha),
      powr_perm = mean(pval_perm[1:rejectH0] < alpha),
      p_sw = rejectH0 / Nsim
    )
  }
  close_cluster(my_cl)
})

# ========================
# Initialize storage arrays
power_t <- power_perm <- prob_sw <- array(NA, dim = c(length(sample_sizes), length(distribution)), 
                               dimnames = list(sample_sizes, distribution))

# Populate results from simulation output
for (i in 1:nrow(params)) {
  n <- as.character(params$n[i])
  dist <- as.character(params$dist[i])
  power_t[n, dist] <- sim_out[[i]]$powr_t
  power_perm[n, dist] <- sim_out[[i]]$powr_perm
  prob_sw[n, dist] <- sim_out[[i]]$p_sw
}

# calculate power loss
powerloss<- power_perm - power_t
# ========================
# Print results
print("Power of t-test:")
print(power_t)

print("Power of permutation test:")
print(power_perm)

print("Power loss for using t-test")
print(powerloss)

expected_powerloss <- round(powerloss * prob_sw, 5)

print("Expected power loss")
print(expected_powerloss)
# ========================
# Optionally save results
save.image(paste0("TwoSample_conditional_powerloss_04152025.RData"))
