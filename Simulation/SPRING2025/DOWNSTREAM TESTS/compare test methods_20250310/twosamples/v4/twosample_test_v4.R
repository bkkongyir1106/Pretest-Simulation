# Set working directory (modify as needed)
base_dir <- "~/Desktop/OSU/Research/Pretest-Simulation"
setwd(file.path(base_dir, "Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/twosamples/v3"))

# Load user-defined functions
source(file.path(base_dir, "functions/utility.R"))
source(file.path(base_dir, "functions/User_defined_functions.R"))


# Simulation parameters
Nsim <- 1e2
sample_size <- c(8, 10, 15, 20, 25, 30, 40, 50)
distribution <- c("Normal", "Exponential", "LogNormal")
alpha <- 0.05
n_boot <- 1e2
effect_size <- 0.0
n_cores <- max(1, detectCores() - 1)

# Initialize result arrays
init_result_array <- function() {
  array(NA, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
}
power.sw.test <- error.t.test <- error.wilcox.test <- error.perm.test <- 
  error.boot.test <- error.conditional.test <- error.t.perm.test <- 
  error.t_perm_split.test <- init_result_array()

# Function to run one simulation iteration
run_simulation <- function(n, dist, alpha, effect_size, n_boot) {
  x <- generate_data(n, dist)
  y <- generate_data(n, dist)
  
  sw_pval <- shapiro.test(x)$p.value
  pval_t <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
  pval_wilcox <- TwoSample.test(x, y, test = "Wilcox", alpha = alpha, effect_size = effect_size)
  pval_perm <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
  pval_boot <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
  
  # Conditional t-test
  pval_conditional <- ifelse(sw_pval > alpha & shapiro.test(y)$p.value > alpha, 
                             t.test(x, y, mu = effect_size)$p.value, NA)
  
  # Two-stage t/perm test
  pval_t_perm <- ifelse(sw_pval > alpha & shapiro.test(y)$p.value > alpha, 
                        TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size), 
                        TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot))
  
  return(c(sw_pval, pval_t, pval_wilcox, pval_perm, pval_boot, pval_conditional, pval_t_perm))
}

# Detect OS type
os_type <- .Platform$OS.type

# Track simulation start time
start_time <- Sys.time()
cat("Simulation started at:", format(start_time), "\n\n")

# Set up progress bar
total_steps <- length(distribution) * length(sample_size)
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
step_counter <- 0

# Loop over conditions
for (i in seq_along(distribution)) {
  dist <- distribution[i]
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    
    if (os_type == "unix") {
      sim_results <- mclapply(1:Nsim, mc.cores = n_cores, function(k) {
        run_simulation(n, dist, alpha, effect_size, n_boot)
      })
    } else {
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      sim_results <- foreach(k = 1:Nsim, .combine = rbind, .packages = c("parallel")) %dopar% {
        source(file.path(base_dir, "functions/User_defined_functions.R"))
        source(file.path(base_dir, "functions/utility.R"))
        run_simulation(n, dist, alpha, effect_size, n_boot)
      }
      stopCluster(cl)
    }
    
    sim_matrix <- do.call(rbind, sim_results)
    
    power.sw.test[j, i] <- mean(sim_matrix[, 1] < alpha)
    error.t.test[j, i] <- mean(sim_matrix[, 2] < alpha)
    error.wilcox.test[j, i] <- mean(sim_matrix[, 3] < alpha)
    error.perm.test[j, i] <- mean(sim_matrix[, 4] < alpha)
    error.boot.test[j, i] <- mean(sim_matrix[, 5] < alpha)
    error.conditional.test[j, i] <- mean(sim_matrix[, 6] < alpha, na.rm = TRUE)
    error.t.perm.test[j, i] <- mean(sim_matrix[, 7] < alpha, na.rm = TRUE)
    
    step_counter <- step_counter + 1
    setTxtProgressBar(pb, step_counter)
  }
}
close(pb)

# Track end time
end_time <- Sys.time()
cat("\n\nSimulation completed at:", format(end_time), "\n")
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# Save results
write.csv(power.sw.test, "power_sw_test.csv")
write.csv(error.t.test, "error_t_test.csv")
write.csv(error.wilcox.test, "error_wilcox_test.csv")
write.csv(error.perm.test, "error_perm_test.csv")
write.csv(error.boot.test, "error_bootstrap_test.csv")
write.csv(error.conditional.test, "error_conditional_test.csv")
write.csv(error.t.perm.test, "error_t_perm_test.csv")

# Save RData
save(power.sw.test, error.t.test, error.wilcox.test, error.perm.test,
     error.boot.test, error.conditional.test, error.t.perm.test,
     file = "twosample.test.methods.RData")

# Plot Curves
# --- Combine all error metrics and label PlotType ---
plot_data <- data.frame(
  SampleSize = rep(sample_size, times = length(distribution) * 8),
  Distribution = rep(rep(distribution, each = length(sample_size)), times = 8),
  ErrorRate = c(as.vector(error.t.test),
                as.vector(error.wilcox.test),
                as.vector(error.perm.test),
                as.vector(error.boot.test),
                as.vector(error.conditional.test),
                as.vector(error.t.perm.test),
                as.vector(error.t_perm_split.test),
                as.vector(power.sw.test)),
  Test = rep(c("t-test", "Wilcox", "Permutation", "Bootstrap",
               "Conditional t-test", "adaptive t/perm", "adaptive t/perm + split", "Shapiro-Wilk"),
             each = length(sample_size) * length(distribution))
)


# Define PlotType for row-based faceting
plot_data$PlotType <- ifelse(
  plot_data$Test == "Shapiro-Wilk", "SW Power",
  ifelse(plot_data$Test == "Conditional t-test", "Conditional Type I Error", "Type I Error")
)

# Convert PlotType to factor for ordering rows
plot_data$PlotType <- factor(plot_data$PlotType, levels = c("SW Power", "Type I Error", "Conditional Type I Error"))

# Base plot
p <- ggplot(plot_data, aes(x = SampleSize, y = ErrorRate, color = Test, linetype = Test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(PlotType ~ Distribution, scales = "free_y") +
  labs(
    title = "Type I Error Rates, Shapiro-Wilk Power, and Conditional Type I Error by Distribution",
    x = "Sample Size",
    y = "Rate / Power"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Add 0.05 reference lines only to Type I Error plots
p <- p + geom_hline(
  data = subset(plot_data, PlotType == "Type I Error"),
  aes(yintercept = 0.05),
  linetype = "dashed",
  color = "gray40"
)

p
# Save plot
ggsave("type1_error_sw_power_conditional_plot.png", p, width = 12, height = 10)
