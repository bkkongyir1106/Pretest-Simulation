# Load user-defined functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/twosamples")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
{
  par_set <- function(cores_reserve = 2) {
    cores <- parallel::detectCores()
    cores_use <- cores - cores_reserve
    if (Sys.info()["sysname"] == "Windows") {
      cl <- parallel::makeCluster(cores_use) # make a socket cluster
      doParallel::registerDoParallel(cl)     # for both Windows & Unix-like
    } else {
      cl <- snow::makeSOCKcluster(cores_use) # make a socket cluster
      doSNOW::registerDoSNOW(cl)             # for non-Windows
    }
    foreach::getDoParWorkers()
    return(cl)
  }
  close_cluster <- function(cl) {
    parallel::stopCluster(cl) # close the cluster
  }
}

{
  Nsim <- 1e3
  sample_size <- c(8, 10, 15, 20, 25, 30)
  distribution <- c("Normal", "Exponential", "LogNormal")
  alpha <- 0.05
  n_boot <- 1e3
  effect_size <- 0.0  
}

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size) * length(distribution) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}

# Run simulation
system.time({
  sim_out <- foreach(i = 1:ntasks,
                     .packages = c("LaplacesDemon", "VGAM", "snow"),
                     .options.snow = opts) %dopar% {
                       n <- sample_size[((i - 1) %% length(sample_size)) + 1]
                       dist <- distribution[((i - 1) %/% length(sample_size)) + 1]
                       set.seed(1234)
                       
      # Storage for p-values for each test method
      pval.t_test <- pval.wilcox.test <- pval.split.test <- pval.perm.test <- pval.boot.test <- pval_t_perm.test <- pval_t_boot.test <- numeric(Nsim)
      for (i in 1:Nsim) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        
        pval.t_test[i] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
        pval.wilcox.test[i] <- TwoSample.test(x, y, test = "Wilcox", alpha = alpha, effect_size = effect_size)
        pval.split.test[i] <- TwoSample.test(x, y, test = "2stage.split", alpha = alpha, effect_size = effect_size)
        pval.perm.test[i] <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
        pval.boot.test[i] <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
        
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
      # Correctly assign results
      Results <- list(
        error.t.test = mean(pval.t_test < alpha),
        error.wilcox.test = mean(pval.wilcox.test < alpha),
        error.split.test = mean(pval.split.test < alpha),
        error.perm.test = mean(pval.perm.test < alpha),
        error.boot.test = mean(pval.boot.test < alpha),
        error.t_perm.test = mean(pval_t_perm.test < alpha),
        error.t_boot.test = mean(pval_t_boot.test < alpha)
      )
      return(Results)
    }
})
close_cluster(my_cl)        

# Initialize each TypeIerror matrix separately
powervec <- numeric(length(sample_size) * length(distribution))
TypeIerror.t.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.wilcox.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.split.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.perm.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.boot.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.t_perm.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.t_boot.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))

# Populate results
for (i in seq_along(sample_size)) {
  for (j in seq_along(distribution)) {
    index <- (j - 1) * length(sample_size) + i
    TypeIerror.t.test[i, j] <- sim_out[[index]]$error.t.test
    TypeIerror.t.test[i, j] <- sim_out[[index]]$error.t.test
    TypeIerror.wilcox.test[i, j] <- sim_out[[index]]$error.wilcox.test
    TypeIerror.split.test[i, j] <- sim_out[[index]]$error.split.test
    TypeIerror.perm.test[i, j] <- sim_out[[index]]$error.perm.test
    TypeIerror.boot.test[i, j] <- sim_out[[index]]$error.boot.test
    TypeIerror.t_perm.test[i, j] <- sim_out[[index]]$error.t_perm.test
    TypeIerror.t_boot.test[i, j] <- sim_out[[index]]$error.t_boot.test
  }
}

# Save the results
save(Nsim, sample_size, TypeIerror.t.test, TypeIerror.wilcox.test, TypeIerror.split.test, TypeIerror.perm.test, 
     TypeIerror.boot.test, TypeIerror.t_perm.test, TypeIerror.t_boot.test,
     file = "two_sample.compare_methods.RData")

# ==============================================================================
# ---------------------- Two sample conditional --------------------------------
# ==============================================================================
# Set seed and simulation parameters

# Parallelized simulation setup
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size) * length(distribution) 
  pb <- txtProgressBar(max=ntasks, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
}

system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM", "snow"),
                     .options.snow=opts) %:%
    foreach(dist = distribution) %dopar% {
      set.seed(1234)
      pvals<- numeric(Nsim)
      iter = 0
      while(iter < Nsim) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        
        if(shapiro.test(x)$p.value > alpha && shapiro.test(y)$p.value > alpha) {
          iter <- iter + 1
          pvals[iter] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = 0, B = NULL)
        }
      }
      # Correctly assign results
      Results <- list(
        error.t.test = mean(pvals < alpha)
    
      )
      return(Results)
    }
})
close_cluster(my_cl)        

# Initialize each TypeIerror matrix separately
powervec <- numeric(length(sample_size) * length(distribution))
conditional.TypeIerror <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))

# Populate results
for (i in seq_along(sample_size)) {
  for (j in seq_along(distribution)) {
    conditional.TypeIerror[i, j] <- sim_out[[i]][[j]]$error.t.test
    
  }
}
# Save the results
save(Nsim, sample_size, conditional.TypeIerror, file = "two_sample.conditional.RData")

# =================================== SW test for normality ====================
# Set up simulation parameters
N <- 1e4
alpha <- 0.05
distribution <- c("Normal", "Exponential", "LogNormal")
test <- "SW"
sample_size <- c(8, 10, 15, 20, 25, 30)

# Parallelized simulation setup
my_cl <- par_set(cores_reserve = 2)
ntasks <- length(sample_size) * length(distribution)
pb <- txtProgressBar(max = ntasks, style = 3)

opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# Perform simulation with parallel processing
system.time({
  sim_out <- foreach(n = sample_size, 
                     .packages = c("LaplacesDemon", "VGAM", "moments", "nortest"), 
                     .options.snow = opts) %:%
    foreach(dist = distribution)%dopar% {
      set.seed(12345)  
      pval <- numeric(N)
      
      for (i in seq_len(N)) {
        x <- generate_data(n, dist)
        pval[i] <- generate_tests(x, test)$p.value
      }
      
      powr <- mean(pval < alpha)
      list(powr = powr)
    }
  
  # Close the parallel cluster
  close_cluster(my_cl)
  
  # Store results into a more efficient structure
  power_norm.test <- array(0, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
  
  # Populate the power results
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distribution)) {
        power_norm.test[i, j] <- sim_out[[i]][[j]]$powr
    }
  }
})

# Output the results
print(power_norm.test)

# Save the results to a file
save.image("sw.test.RData")



# ===================== plots =================================================
# Load the datasets
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/onesample/one_sample.compare_methods.RData")
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/onesample/sw.test.RData")
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/onesample/one_sample.conditional.RData")
typeI_data <- data.frame(
  sample_size = sample_size,
  Exp_t = TypeIerror.t.test[, "Exponential"],
  Exp_split = TypeIerror.split.test[, "Exponential"],
  Exp_boot = TypeIerror.boot.test[, "Exponential"],
  Exp_t_boot = TypeIerror.t_boot.test[, "Exponential"],
  LogNorm_t = TypeIerror.t.test[, "LogNormal"],
  LogNorm_split = TypeIerror.split.test[, "LogNormal"],
  LogNorm_boot = TypeIerror.boot.test[, "LogNormal"],
  LogNorm_t_boot = TypeIerror.t_boot.test[, "LogNormal"]
)

# Reshape the data for plotting
typeI_long <- typeI_data %>%
  pivot_longer(
    cols = starts_with("Exp_") | starts_with("LogNorm_"),
    names_to = "test_method",
    values_to = "typeI_error"
  ) %>%
  mutate(
    distribution = ifelse(grepl("^Exp_", test_method), "Exponential", "LogNormal"),
    test_method = gsub("^Exp_|^LogNorm_", "", test_method)
  )

# Combine power data into a data frame
power_data <- data.frame(
  sample_size = sample_size,
  Exponential = power_norm.test[, "Exponential"],
  LogNormal = power_norm.test[, "LogNormal"]
)

# Reshape the power data for plotting
power_long <- power_data %>%
  pivot_longer(
    cols = c(Exponential, LogNormal),
    names_to = "distribution",
    values_to = "power"
  )

# Combine power data into a data frame
conditional_typeI_df <- data.frame(
  sample_size = sample_size,
  Exponential = conditional.TypeIerror[, "Exponential"],
  LogNormal = conditional.TypeIerror[, "LogNormal"]
)

# Reshape the power data for plotting
conditional_typeI_long <- conditional_typeI_df %>%
  pivot_longer(
    cols = c(Exponential, LogNormal),
    names_to = "distribution",
    values_to = "conditional_typeI_error"
  )

# Plot conditional Type I error rates for Exponential and LogNormal
conditional_typeI_plot <- ggplot(conditional_typeI_long, aes(x = sample_size, y = conditional_typeI_error, color = distribution)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylim(0, 0.25) +  # Set y-axis limits for conditional Type I error
  facet_wrap(~distribution, scales = "free_y") +
  labs(
    x = "Sample Size",
    y = "Type I error",
    title = "Conditional Type I error",
    color = "Distribution"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot Type I error rates for all test methods in the same plot for each distribution
typeI_plot <- ggplot(typeI_long, aes(x = sample_size, y = typeI_error, color = test_method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylim(0, 0.07) +
  facet_wrap(~distribution, scales = "free_y") +
  labs(
    x = "Sample Size",
    y = "Type I Error Rate",
    title = "Type I Error Rates by Test Method and Distribution",
    color = "Test Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot power of the Shapiro-Wilk test for each distribution
power_plot <- ggplot(power_long, aes(x = sample_size, y = power, color = distribution)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~distribution, scales = "free_y") +
  labs(
    x = "Sample Size",
    y = "Power of SW Test",
    title = "Power of Shapiro-Wilk Test by Distribution",
    color = "Distribution"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine the plots vertically
combined_plot <- conditional_typeI_plot / typeI_plot / power_plot +
  plot_layout(heights = c(1, 2, 1))  # Adjust heights as needed

# Display the combined plot
print(combined_plot)

