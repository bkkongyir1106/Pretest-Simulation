# Load user-defined functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/twosamples")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
{
  Nsim <- 1e3
  sample_size <- c(8, 10, 15, 20, 25, 30)
  distribution <- c("Normal", "Exponential", "LogNormal")
  alpha <- 0.05
  n_boot <- 1e3
  effect_size <- 0.0  
}

# define a function for arrays  
init_result_array <- function() {
  array(data = NA, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
}
# create arrays for each test method
power.sw.test <- error.t.test <- error.wilcox.test <- error.perm.test <- error.boot.test <- error.t_perm_split.test <- init_result_array()

for(i in seq_along(distribution)){
  dist <- distribution[i]
  for(j in seq_along(sample_size)){
    n <- sample_size[j]
    pval.sw.test <- pval.t_test <- pval.wilcox.test <- pval.split.test <- pval.perm.test <- pval.boot.test <- numeric(Nsim)
    
    for(k in 1:Nsim){
      x <- generate_data(n, dist)
      y <- generate_data(n, dist)
      
      # SW test
      pval.sw.test[k] <- shapiro.test(x)$p.value
      
      # Main tests
      pval.t_test[k] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size)
      pval.wilcox.test[k] <- TwoSample.test(x, y, test = "Wilcox", alpha = alpha, effect_size = effect_size)
      pval.perm.test[k] <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      pval.boot.test[k] <- two_sample_bootstrap_p_value(x, y, effect_size = effect_size, n_bootstrap = n_boot)
      
      # Split-sample test
      split_point <- floor(length(x) / 2)
      first_half.x <- x[1:split_point]
      second_half.x <- x[(split_point + 1):length(x)]
      first_half.y <- y[1:split_point]
      second_half.y <- y[(split_point + 1):length(y)]  
      
      if(shapiro.test(first_half.x)$p.value > alpha & shapiro.test(first_half.y)$p.value > alpha){
        pval.split.test[k] <- t.test(second_half.x, second_half.y, mu = effect_size)$p.value
      } else {
        pval.split.test[k] <- TwoSample.test(second_half.x, second_half.y, test = "perm", alpha = alpha, effect_size = effect_size, B = n_boot)
      }
    }
    
    # Save results
    power.sw.test[j, i] <- mean(pval.sw.test < alpha)
    error.t.test[j, i] <- mean(pval.t_test < alpha)
    error.wilcox.test[j, i] <- mean(pval.wilcox.test < alpha)
    error.perm.test[j, i] <- mean(pval.perm.test < alpha)
    error.boot.test[j, i] <- mean(pval.boot.test < alpha)
    error.t_perm_split.test[j, i] <- mean(pval.split.test < alpha)
  }
}

# Save results to CSV
write.csv(power.sw.test, "power_sw_test.csv")
write.csv(error.t.test, "error_t_test.csv")
write.csv(error.wilcox.test, "error_wilcox_test.csv")
write.csv(error.perm.test, "error_perm_test.csv")
write.csv(error.boot.test, "error_bootstrap_test.csv")
write.csv(error.t_perm_split.test, "error_t_perm_split_test.csv")

# Combine all error rate matrices into a single data frame for plotting
plot_data <- data.frame(
  SampleSize = rep(sample_size, times = length(distribution) * 6),
  Distribution = rep(rep(distribution, each = length(sample_size)), times = 6),
  ErrorRate = c(as.vector(error.t.test),
                as.vector(error.wilcox.test),
                as.vector(error.perm.test),
                as.vector(error.boot.test),
                as.vector(error.t_perm_split.test),
                as.vector(power.sw.test)),
  Test = rep(c("t-test", "Wilcox", "Permutation", "Bootstrap", "T+Split", "Shapiro-Wilk"), each = length(sample_size) * length(distribution))
)

# Plot
ggplot(plot_data, aes(x = SampleSize, y = ErrorRate, color = Test, linetype = Test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Distribution, scales = "free_y") +
  labs(
    title = "Type I Error Rates / Power vs. Sample Size by Distribution",
    x = "Sample Size",
    y = "Error Rate / SW Power"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Save plot as image
ggsave("type1_error_rates_plot.png", width = 10, height = 6)
