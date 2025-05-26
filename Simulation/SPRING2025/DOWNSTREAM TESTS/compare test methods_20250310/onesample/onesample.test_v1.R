# Load necessary libraries
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)

# Load user-defined functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/onesample")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Simulation parameters
Nsim <- 1e4
sample_size <- c(8, 10, 15, 20, 25, 30, 40, 50)
distribution <- c("Normal", "Exponential", "LogNormal")
alpha <- 0.05
n_boot <- 1e4
effect_size <- 0.0
n_cores <- max(1, detectCores() - 1)

# Initialize result arrays
init_result_array <- function() {
  array(data = NA, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
}
power.sw.test <- error.t.test <- error.wilcox.test <- error.boot.test <- error.conditional.test <- error.t.boot.test <- error.t_boot_split.test <- init_result_array()

# Function to run one simulation iteration
run_simulation <- function(n, dist, alpha, effect_size, n_boot) {
  x <- generate_data(n, dist)
 
  sw_pval <- shapiro.test(x)$p.value
  pval_t <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size)
  pval_wilcox <- OneSample.test(x, test = "Wilcox", alpha = alpha, effect_size = effect_size)
  pval_boot <- one_sample_bootstrap_p_value(x, effect_size = effect_size, n_bootstrap = n_boot)
  # conditional t-test
  if (shapiro.test(x)$p.value > alpha) {
    pval_conditional <- t.test(x, mu = effect_size)$p.value
  } else {
    pval_conditional <- NA
  }
  
  # two-stage t/perm test
  if (shapiro.test(x)$p.value > alpha ) {
    pval_t_boot <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size, B = NULL)
  } else {
    pval_t_boot <- one_sample_bootstrap_p_value(x, effect_size = effect_size, n_bootstrap = n_boot)
  }
  # conditional split t/perm test
  split_point <- floor(length(x) / 2)
  first_half.x <- x[1:split_point]
  second_half.x <- x[(split_point + 1):length(x)]
  
  if (shapiro.test(first_half.x)$p.value > alpha) {
    pval_split <- t.test(second_half.x, mu = effect_size)$p.value
  } else {
    pval_split <- one_sample_bootstrap_p_value(second_half.x, effect_size = effect_size, n_bootstrap = n_boot)
  }
  
  return(c(sw_pval, pval_t, pval_wilcox, pval_boot, pval_conditional, pval_t_boot,  pval_split))
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
      # macOS/Linux: use mclapply
      sim_results <- mclapply(1:Nsim, mc.cores = n_cores, function(k) {
        run_simulation(n, dist, alpha, effect_size, n_boot)
      })
    } else {
      # Windows: use foreach + doParallel
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      sim_results <- foreach(k = 1:Nsim, .combine = rbind, .packages = c()) %dopar% {
        source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
        source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
        run_simulation(n, dist, alpha, effect_size, n_boot)
      }
      stopCluster(cl)
    }
    
    sim_matrix <- do.call(rbind, sim_results)
    
    power.sw.test[j, i] <- mean(sim_matrix[, 1] < alpha)
    error.t.test[j, i] <- mean(sim_matrix[, 2] < alpha)
    error.wilcox.test[j, i] <- mean(sim_matrix[, 3] < alpha)
    error.boot.test[j, i] <- mean(sim_matrix[, 4] < alpha)
    error.conditional.test[j, i] <- round(mean(sim_matrix[, 5] < alpha, na.rm = TRUE), 3)
    error.t.boot.test[j, i] <- mean(sim_matrix[, 6] < alpha, na.rm = TRUE)
    error.t_boot_split.test[j, i] <- mean(sim_matrix[, 7] < alpha)
    # Update progress bar
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
write.csv(error.boot.test, "error_bootstrap_test.csv")
write.csv(error.t_boot_split.test, "error_t_boot_split_test.csv")
write.csv(error.conditional.test, "error_conditional_test.csv")
write.csv(error.t.boot.test, "error_t_boot_test.csv")

# save RData
save(power.sw.test, error.t.test, error.wilcox.test, error.boot.test, error.t_boot_split.test,error.conditional.test, error.t_boot_split.test, file = "twosample.test.methods.RData" )
# Plot Curves
# --- Combine all error metrics and label PlotType ---
plot_data <- data.frame(
  SampleSize = rep(sample_size, times = length(distribution) * 7),
  Distribution = rep(rep(distribution, each = length(sample_size)), times = 7),
  ErrorRate = c(as.vector(error.t.test),
                as.vector(error.wilcox.test),
                as.vector(error.boot.test),
                as.vector(error.conditional.test),
                as.vector(error.t.boot.test),
                as.vector(error.t_boot_split.test),
                as.vector(power.sw.test)),
  Test = rep(c("t-test", "Wilcox", "Bootstrap",
               "Conditional t-test", "t/boot", "t/boot+Split", "Shapiro-Wilk"),
             each = length(sample_size) * length(distribution))
)

# Filter out Normal distribution
#plot_data <- subset(plot_data, Distribution != "Normal")

# Define PlotType for row-based faceting
plot_data$PlotType <- ifelse(
  plot_data$Test == "Shapiro-Wilk", "SW Power",
  ifelse(plot_data$Test == "Conditional t-test", "Conditional Type I Error", "Type I Error")
)

# Convert PlotType to factor for ordering rows
plot_data$PlotType <- factor(plot_data$PlotType, levels = c("SW Power", "Type I Error", "Conditional Type I Error"))

# --- Plot using facet_grid (rows = PlotType, columns = Distribution) ---
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

# Add 0.05 reference lines **only to Type I Error plots**
# We'll layer geom_hline() for those panels using data-specific filtering
p <- p + geom_hline(
  data = subset(plot_data, PlotType == "Type I Error"),
  aes(yintercept = 0.05),
  linetype = "dashed",
  color = "gray40"
)

p
# Save plot
ggsave("type1_error_sw_power_conditional_plot.png", p, width = 12, height = 10)
