# Clear environment
rm(list = ls())

# Set working directory and source custom functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE/TYPE I ERROR")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

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

# Parameters
N <- 1e4          # Number of simulations
alpha <- 0.05     # Significance level
distribution <- c("Normal", "Exponential")
sample_sizes <- c(8, 10, 15, 20, 25, 30, 50)

# Setup parallel backend and progress bar
my_cl <- par_set(cores_reserve = 2)
params <- expand.grid(n = sample_sizes, dist = distribution)
ntasks <- nrow(params)
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run simulation
system.time({
  sim_out <- foreach(i = 1:nrow(params), .packages = c("LaplacesDemon", "VGAM"), .options.snow = opts) %dopar% {
    n <- params$n[i]
    dist <- params$dist[i]
    
    set.seed(12345 + i)  # Ensure different seeds per job
    Nsim <- 0
    sw_not_rejected <- 0
    pval <- numeric(N)
    while ( sw_not_rejected < N) {
      x <- generate_data(n, dist)

      Nsim <- Nsim + 1
      
      # Conditional adaptive test
      if (shapiro.test(x)$p.value > alpha) {
        sw_not_rejected <- sw_not_rejected + 1
        pval[sw_not_rejected] <- t.test(x)$p.value
      } 
    }
    
    list(error = mean(pval < alpha))
    
    
  }
  close_cluster(my_cl)
})

# Initialize storage arrays
conditionalTypeI.errorRate <- array(NA, dim = c(length(sample_sizes), length(distribution)), dimnames = list(sample_sizes, distribution))

# Populate results from simulation output
for (i in 1:nrow(params)) {
  n <- as.character(params$n[i])
  dist <- as.character(params$dist[i])
  conditionalTypeI.errorRate[n, dist] <- sim_out[[i]]$error
}

# print
print(conditionalTypeI.errorRate)

# Save data
save.image(paste0("OneSampleconditionalTypeI_errorRate_04132025",".RData"))
