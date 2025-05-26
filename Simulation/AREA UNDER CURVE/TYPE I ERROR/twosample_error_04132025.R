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
N <- 1e4          
B <- 1e3          
alpha <- 0.05     
distribution <- c("Normal", "Chi-Square", "LogNormal", "Exponential")
sample_sizes <- c(8, 10, 15, 20, 25, 30, 50)

# Test statistic: 
calculate_test_statistic <- function(x, y) {
  (mean(x) - mean(y)) / sqrt(var(x)/length(x) + var(y)/length(y))
}

# Trapezoidal integration for area under curve
compute_area <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(sample_sizes) - min(sample_sizes))
}

# Setup parallel and progress bar
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
    
    set.seed(12345 + i)  
    pval_t <- pval_wilcox <- pval_t_or_wilcox <- pval_t_or_perm<- pval_perm <- numeric(N)
    pval_sw_x <- pval_sw_y <- numeric(N)
    for (s in 1:N) {
      x <- generate_data(n, dist)
      y <- generate_data(n, dist)
      
      # shapiro-wilk test
      pval_sw_x[s] <- shapiro.test(x)
      pval_sw_y[s] <- shapiro.test(y)
      
      # Classical tests
      pval_t[s] <- t.test(x, y)$p.value
      pval_wilcox[s] <- wilcox.test(x, y)$p.value
      
      # Permutation test
      data <- c(x, y)
      observe_stat <- calculate_test_statistic(x, y)
      permuted_stat <- numeric(B)
      
      for (j in 1:B) {
        permuted <- sample(data)
        sample_x <- permuted[1:n]
        sample_y <- permuted[(n + 1):(2 * n)]
        permuted_stat[j] <- calculate_test_statistic(sample_x, sample_y)
      }
      pval_perm[s] <- mean(abs(permuted_stat) >= abs(observe_stat))
    }
  
    list(
      error_t = mean(pval_t < alpha),
      error_wilcox = mean(pval_wilcox < alpha),
      error_perm = mean(pval_perm < alpha)
    )
  }
  close_cluster(my_cl)
})

# Initialize storage arrays
dims <- c(length(sample_sizes), length(distribution))
dimnames_list <- list(sample_sizes, distribution)
TypeI.errorRate_t <- TypeI.errorRate_wilcox <- TypeI.errorRate_t_wilcox <- TypeI.errorRate_t_perm<- TypeI.errorRate_perm <-
  array(NA, dim = dims, dimnames = dimnames_list)

# Populate results 
for (i in 1:nrow(params)) {
  n <- as.character(params$n[i])
  dist <- as.character(params$dist[i])
  TypeI.errorRate_t[n, dist]         <- sim_out[[i]]$error_t
  TypeI.errorRate_wilcox[n, dist]    <- sim_out[[i]]$error_wilcox
  TypeI.errorRate_perm[n, dist]      <- sim_out[[i]]$error_perm
}

# Compute AUC for each test
area_t        <- apply(TypeI.errorRate_t, 2, compute_area, x = sample_sizes)
area_wilcox   <- apply(TypeI.errorRate_wilcox, 2, compute_area, x = sample_sizes)
area_perm     <- apply(TypeI.errorRate_perm, 2, compute_area, x = sample_sizes)

# Output
print(TypeI.errorRate_t)
print(TypeI.errorRate_wilcox)
print(TypeI.errorRate_perm)

print(area_t)
print(area_wilcox)
print(area_perm)

# Save data
#save.image(paste0("TwoSampleAUCTypeI_errorRate_04132025",".RData"))
