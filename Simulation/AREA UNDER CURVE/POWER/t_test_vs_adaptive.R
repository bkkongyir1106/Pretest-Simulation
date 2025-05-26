# Set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/AREA UNDER CURVE/POWER")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

{
  Nsim <- 1e4
  B <- 1e3
  alpha <- 0.05
  distributions <- c("Normal", "Uniform","Exponential",  "LogNormal")
  #distributions <- c("Normal", "LogNormal")
  sample_sizes <- c(8, 10, 15, 20, 25, 30, 50)
  effect_size <- 0.5
}

calculate_test_statistic <- function(x, y) {
  (mean(x) - mean(y)) / sqrt(var(x)/length(x) + var(y)/length(y))
}

permutation_test <- function(x, y, B = 1000) {
  observed <- calculate_test_statistic(x, y)
  combined <- c(x, y)
  perm_stats <- replicate(B, {
    permuted <- sample(combined)
    x_star <- permuted[1:length(x)]
    y_star <- permuted[(length(x)+1):(2*length(x))]
    calculate_test_statistic(x_star, y_star)
  })
  mean(abs(perm_stats) >= abs(observed))
}

run_simulation <- function(dist, n, Nsim, effect_size, alpha, B) {
  reject_perm_null <- reject_perm_alt <- numeric(Nsim)
  reject_adapt_null <- reject_adapt_alt <- numeric(Nsim)
  reject_t_null <- reject_t_alt <- numeric(Nsim)
  N_t = 0
  N_sw = 0
  for (i in 1:Nsim) {
    # Under null
    x0 <- generate_data(n, dist)
    y0 <- generate_data(n, dist)
    N_t = N_t + 1
    pval_t_null <- t.test(x0, y0)$p.value
    reject_t_null[i] <- (pval_t_null < alpha)
    
    p_null <- permutation_test(x0, y0, B)
    reject_perm_null[i] <- (p_null < alpha)
    
    if (shapiro.test(x0)$p.value > alpha && shapiro.test(y0)$p.value > alpha) {
      N_sw = N_sw + 1
      p_t_null <- t.test(x0, y0)$p.value
      reject_adapt_null[i] <- (p_t_null < alpha)
    } else {
      reject_adapt_null[i] <- (p_null < alpha)
    }
    
    # Under alternative
    x1 <- x0
    y1 <- y0 + effect_size
    
    pval_t_alt <- t.test(x0, y0 + effect_size)$p.value
    reject_t_alt[i] <- (pval_t_alt < alpha)
    
    p_alt <- permutation_test(x1, y1, B)
    reject_perm_alt[i] <- (p_alt < alpha)
    
    if (shapiro.test(x1)$p.value > alpha && shapiro.test(y1)$p.value > alpha) {
      p_t_alt <- t.test(x1, y1)$p.value
      reject_adapt_alt[i] <- (p_t_alt < alpha)
    } else {
      reject_adapt_alt[i] <- (p_alt < alpha)
    }
  }
  
  list(
    type1_perm = mean(reject_perm_null),
    type1_t = mean(reject_t_null),
    type1_adapt = mean(reject_adapt_null),
    power_t = mean(reject_t_alt),
    power_perm = mean(reject_perm_alt),
    power_adapt = mean(reject_adapt_alt),
    p_sw = N_sw/N_t
  )
}

# ------------ Initialize Results -----------#
dims <- c(length(sample_sizes), length(distributions))
dimnames_list <- list(as.character(sample_sizes), distributions)

power_t_mat <- power_perm_mat <- power_adapt_mat <- array(NA, dim = dims, dimnames = dimnames_list)
type1_t_mat <- type1_perm_mat <- type1_adapt_mat <- prob_sw <- array(NA, dim = dims, dimnames = dimnames_list)

# ---------- Parallel Setup ---------------- #
cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(sample_sizes) * length(distributions), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

params <- expand.grid(n = sample_sizes, dist = distributions)

results <- foreach(i = 1:nrow(params), .options.snow = opts, .packages = c()) %dopar% {
  n <- params$n[i]
  dist <- params$dist[i]
  run_simulation(dist, n, Nsim, effect_size, alpha, B)
}

close(pb)
stopCluster(cl)

# ----------- Store Results -------- #
for (i in 1:nrow(params)) {
  n <- as.character(params$n[i])
  dist <- as.character(params$dist[i])
  res <- results[[i]]
  
  power_t_mat[n, dist] <- res$power_t
  power_perm_mat[n, dist] <- res$power_perm
  power_adapt_mat[n, dist] <- res$power_adapt
  type1_t_mat[n, dist] <- res$type1_t
  type1_perm_mat[n, dist] <- res$type1_perm
  type1_adapt_mat[n, dist] <- res$type1_adapt
  prob_sw[n, dist] <- res$p_sw
}

power_diff <- power_perm_mat - power_t_mat
powerloss <- power_adapt_mat - power_t_mat
inflation_error <- type1_t_mat - alpha

Expected_powerloss <- round(powerloss * (1 - prob_sw), 5)
Expectedinflation_error <- round(inflation_error *  (1 - prob_sw), 5)
# ----------- Print Results ----------- #
cat("\nPower (t-test):\n"); print(power_t_mat)
cat("\nPower (Permutation):\n"); print(power_perm_mat)
cat("\nPower (Adaptive):\n"); print(power_adapt_mat)
cat("\nPower (power_diff):\n"); print(power_diff)

cat("\nType I Error (Permutation):\n"); print(type1_perm_mat)
cat("\nType I Error (t-test):\n"); print(type1_t_mat)
cat("\nType I Error (Adaptive):\n"); print(type1_adapt_mat)
cat("\nType I Error (inflation_error):\n"); print(inflation_error)
cat("\nType I Error (Expected_inflation):\n"); print(Expectedinflation_error)
cat("\nprob_sw:\n"); print(prob_sw)

# ------- Save for Later Use ----------- #
#save.image("t_test_vs_Adaptive_Simulation.RData")
save.image("t_test_vs_Adaptive_Simulation_pA.RData")

