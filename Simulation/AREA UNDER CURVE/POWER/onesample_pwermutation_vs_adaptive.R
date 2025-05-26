# Set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE/POWER")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Parameters
Nsim <- 1000
B <- 1000
alpha <- 0.05
distributions <- c("Normal", "Uniform", "Exponential", "LogNormal")
sample_sizes <- c(8, 10, 15, 20, 25, 30, 50)
effect_size <- 0.5

# One-sample test statistic (testing mean = 0)
calculate_test_statistic <- function(x, mu0 = 0) {
  (mean(x) - mu0) / (sd(x) / sqrt(length(x)))
}

# One-sample permutation test
permutation_test <- function(x, B = 1000, mu0 = 0) {
  observed <- calculate_test_statistic(x, mu0)
  perm_stats <- replicate(B, {
    signs <- sample(c(-1, 1), length(x), replace = TRUE)
    calculate_test_statistic(signs * x, mu0)
  })
  mean(abs(perm_stats) >= abs(observed))
}

# Simulation function
run_simulation <- function(dist, n, Nsim, effect_size, alpha, B) {
  reject_perm_null <- reject_perm_alt <- numeric(Nsim)
  reject_adapt_null <- reject_adapt_alt <- numeric(Nsim)
  N_t <- N_sw <- 0
  
  for (i in 1:Nsim) {
    # Under null
    x0 <- generate_data(n, dist)
    N_t <- N_t + 1
    p_null <- permutation_test(x0, B)
    reject_perm_null[i] <- (p_null < alpha)
    
    if (shapiro.test(x0)$p.value > alpha) {
      N_sw <- N_sw + 1
      p_t_null <- t.test(x0, mu = 0)$p.value
      reject_adapt_null[i] <- (p_t_null < alpha)
    } else {
      reject_adapt_null[i] <- (p_null < alpha)
    }
    
    # Under alternative
    x1 <- x0 + effect_size
    p_alt <- permutation_test(x1, B)
    reject_perm_alt[i] <- (p_alt < alpha)
    
    if (shapiro.test(x0)$p.value > alpha) {
      p_t_alt <- t.test(x1, mu = 0)$p.value
      reject_adapt_alt[i] <- (p_t_alt < alpha)
    } else {
      reject_adapt_alt[i] <- (p_alt < alpha)
    }
  }
  
  list(
    type1_perm = mean(reject_perm_null),
    type1_adapt = mean(reject_adapt_null),
    power_perm = mean(reject_perm_alt),
    power_adapt = mean(reject_adapt_alt),
    p_sw = N_sw / N_t
  )
}

# Result arrays
dims <- c(length(sample_sizes), length(distributions))
dimnames_list <- list(as.character(sample_sizes), distributions)

power_perm_mat <- power_adapt_mat <- array(NA, dim = dims, dimnames = dimnames_list)
type1_perm_mat <- type1_adapt_mat <- prob_sw <- array(NA, dim = dims, dimnames = dimnames_list)

# Parallel setup
library(doSNOW)
library(foreach)
cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(sample_sizes) * length(distributions), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

params <- expand.grid(n = sample_sizes, dist = distributions)

results <- foreach(i = 1:nrow(params), .options.snow = opts) %dopar% {
  n <- params$n[i]
  dist <- params$dist[i]
  run_simulation(dist, n, Nsim, effect_size, alpha, B)
}

close(pb)
stopCluster(cl)

# Store results
for (i in 1:nrow(params)) {
  n <- as.character(params$n[i])
  dist <- as.character(params$dist[i])
  res <- results[[i]]
  
  power_perm_mat[n, dist] <- res$power_perm
  power_adapt_mat[n, dist] <- res$power_adapt
  type1_perm_mat[n, dist] <- res$type1_perm
  type1_adapt_mat[n, dist] <- res$type1_adapt
  prob_sw[n, dist] <- res$p_sw
}

powerloss <- power_perm_mat - power_adapt_mat
Expected_powerloss <- powerloss * prob_sw

# Print results
cat("\nPower (Permutation):\n"); print(power_perm_mat)
cat("\nPower (Adaptive):\n"); print(power_adapt_mat)
cat("\nType I Error (Permutation):\n"); print(type1_perm_mat)
cat("\nType I Error (Adaptive):\n"); print(type1_adapt_mat)
cat("\nProportion Normal (SW Passed):\n"); print(prob_sw)

# Optionally save
# save.image("OneSample_Permutation_vs_Adaptive_Simulation.RData")
