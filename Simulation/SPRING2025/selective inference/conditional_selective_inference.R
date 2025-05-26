# Settings
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
set.seed(12345)
mu0 <- 0
alpha <- 0.05
alpha_pre <- 0.05
M <- 100000  # Number of replications
sample_size <- c(8, 10, 15, 20, 25, 30)
distribution <- c("Normal", "LogNormal", "Exponential")

# Initialize counters
reject_unconditional <- 0
reject_conditional_naive <- 0
reject_conditional_selective <- 0
sw_pass_count <- 0

# Storage for simulation
t_vals_all <- numeric(M)
sw_pvals_all <- numeric(M)
for(j in seq_along(sample_size)){
  n <- sample_size[j]
  for(k in seq_along(distribution)){
    dist <- distribution[k]
    
    # Step 1: Simulate M datasets under H0:
    for (i in 1:M) {
      x <- generate_data(n, dist)
      t_stat <- (mean(x) - mu0) / (sd(x) / sqrt(n))
      sw_pval <- shapiro.test(x)$p.value
      
      # Unconditional t-test
      if (abs(t_stat) > qt(1 - alpha / 2, df = n - 1)) {
        reject_unconditional <- reject_unconditional + 1
      }
      
      # Conditional t-test (naive, without adjustment)
      if (sw_pval > alpha_pre) {
        sw_pass_count <- sw_pass_count + 1
        if (abs(t_stat) > qt(1 - alpha / 2, df = n - 1)) {
          reject_conditional_naive <- reject_conditional_naive + 1
        }
      }
      
      # Store for conditional selective inference
      t_vals_all[i] <- t_stat
      sw_pvals_all[i] <- sw_pval
    }
    
    # Step 2: Estimate conditional null distribution
    t_cond <- t_vals_all[sw_pvals_all > alpha_pre]
    crit_val_selective <- quantile(abs(t_cond), probs = 1 - alpha / 2)
    
    # Step 3: Apply conditional test with corrected critical value
    for (i in which(sw_pvals_all > alpha_pre)) {
      if (abs(t_vals_all[i]) > crit_val_selective) {
        reject_conditional_selective <- reject_conditional_selective + 1
      }
    }
    
  }
}
# Step 4: Calculate error rates
type1_unconditional <- reject_unconditional / M
type1_conditional_naive <- reject_conditional_naive / sw_pass_count
type1_conditional_selective <- reject_conditional_selective / sw_pass_count
sw_pass_rate <- sw_pass_count / M

# Results
cat("Results from", M, "replications with n =", n, "\n")
cat("-------------------------------------------------\n")
cat("Unconditional Type I Error (no pretest):", round(type1_unconditional, 4), "\n")
cat("Conditional Type I Error (naive t-test):", round(type1_conditional_naive, 4), "\n")
cat("Conditional Type I Error (selective):    ", round(type1_conditional_selective, 4), "\n")
cat("Rate of passing Shapiro–Wilk test:", round(sw_pass_rate, 4), "\n")
cat("Selective critical value (|T|):", round(crit_val_selective, 4), "\n")
