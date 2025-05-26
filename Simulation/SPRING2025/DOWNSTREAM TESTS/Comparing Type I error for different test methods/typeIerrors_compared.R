# Load user-defined functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/Comparing Type I error for different test methods")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

# Set seed and simulation parameters
set.seed(12345)
Nsim <- 100
sample_size_vec <- c(8, 10)#, 15, 20, 25, 30)
distribution_vec <- c("Normal", "Exponential", "LogNormal")
alpha <- 0.05
n_boot <- 100
effect_size <- 0.0  # Type I error = Null is true

# Initialize result matrices
t_test_mat     <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))
wilcox_mat     <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))
split_mat      <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))
perm_mat       <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))
boot_mat       <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))
t_perm_mat     <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))
t_boot_mat     <- matrix(NA, nrow = length(sample_size_vec), ncol = length(distribution_vec))

rownames(t_test_mat) <- sample_size_vec
colnames(t_test_mat) <- distribution_vec
rownames(wilcox_mat) <- rownames(split_mat) <- rownames(perm_mat) <- rownames(boot_mat) <- rownames(t_perm_mat) <- rownames(t_boot_mat) <- sample_size_vec
colnames(wilcox_mat) <- colnames(split_mat) <- colnames(perm_mat) <- colnames(boot_mat) <- colnames(t_perm_mat) <- colnames(t_boot_mat) <- distribution_vec

handlers(global = TRUE)
handlers("txtprogressbar")  # Ensure progress bar works in non-parallel runs

with_progress({
  p <- progressor(steps = length(distribution_vec) * length(sample_size_vec))
# Start simulation
for (d in seq_along(distribution_vec)) {
  dist_current <- distribution_vec[d]
  
  for (s in seq_along(sample_size_vec)) {
    n <- sample_size_vec[s]
    
    # Storage for p-values for each test method
    pval.t_test <- numeric(Nsim)
    pval.wilcox.test <- numeric(Nsim)
    pval.split.test <- numeric(Nsim)
    pval.perm.test <- numeric(Nsim)
    pval.boot.test <- numeric(Nsim)
    pval_t_perm.test <- numeric(Nsim)
    pval_t_boot.test <- numeric(Nsim)
    
    for (i in 1:Nsim) {
      x <- generate_data(n, dist_current)
      y <- generate_data(n, dist_current)
      
      pval.t_test[i]     <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      pval.wilcox.test[i]<- TwoSample.test(x, y, test = "Wilcox", alpha = alpha, effect_size = effect_size)
      pval.split.test[i] <- TwoSample.test(x, y, test = "2stage.split", alpha = alpha, effect_size = effect_size)
      pval.perm.test[i]  <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      pval.boot.test[i]  <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
      
      # Adaptive t/permutation
      if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
        pval_t_perm.test[i] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      } else {
        pval_t_perm.test[i] <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      }
      
      # Adaptive t/bootstrap
      if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
        pval_t_boot.test[i] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      } else {
        pval_t_boot.test[i] <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
      }
    }
    
    # Fill in the type I error matrices
    t_test_mat[s, d]   <- mean(pval.t_test < alpha)
    wilcox_mat[s, d]   <- mean(pval.wilcox.test < alpha)
    split_mat[s, d]    <- mean(pval.split.test < alpha)
    perm_mat[s, d]     <- mean(pval.perm.test < alpha)
    boot_mat[s, d]     <- mean(pval.boot.test < alpha)
    t_perm_mat[s, d]   <- mean(pval_t_perm.test < alpha)
    t_boot_mat[s, d]   <- mean(pval_t_boot.test < alpha)
    
    # Simulated computations (Replace with actual test code)
    Sys.sleep(0.1)  # Simulates computation time
    p()  # Increment progress bar
  }
}
})

# Print the results nicely
cat("\nType I Error - t-test:\n")
print(round(t_test_mat, 3))

cat("\nType I Error - Wilcoxon test:\n")
print(round(wilcox_mat, 3))

cat("\nType I Error - 2-stage split test:\n")
print(round(split_mat, 3))

cat("\nType I Error - Permutation test:\n")
print(round(perm_mat, 3))

cat("\nType I Error - Bootstrap test:\n")
print(round(boot_mat, 3))

cat("\nType I Error - Adaptive t/perm test:\n")
print(round(t_perm_mat, 3))

cat("\nType I Error - Adaptive t/bootstrap test:\n")
print(round(t_boot_mat, 3))

