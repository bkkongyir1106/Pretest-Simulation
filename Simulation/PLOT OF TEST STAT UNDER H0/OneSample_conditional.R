rm(list = ls())
library(ggplot2)
library(gridExtra)
setwd("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/PLOT OF TEST STAT UNDER H0")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

dist_sum <- c("Standard Normal", "Uniform", "t", "Exponential", "Chi-Square", "LogNormal")
set.seed(12345)
N <- 1e5
alpha <- 0.05
d <- 0.5
n <- 20
df <- n - 1

# Initialize list for storing statistics
test_statistic <- Type1_error <- vector("list", length(dist_sum))

for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  Nsim <- 0
  Npassed <- 0
  stat <- numeric(N)  # Preallocate vector for statistics
  
  while (Npassed < N) {
    x <- generate_data(n, dist)
    Nsim <- Nsim + 1
    if (shapiro.test(x)$p.value > alpha) {
      Npassed <- Npassed + 1
      stat[Npassed] <- t.test(x)$statistic
    }
  }
  test_statistic[[i]] <- stat
  #normal Distn
  left_norm <- mean(test_statistic[[i]] < qt(0.025, df = df))
  right_norm <- mean(test_statistic[[i]] > qt(0.975, df = df))
  Type1_error[[i]] <- sum(left_norm + right_norm)
  
}

# Output each statistic by distribution name
for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  cat("Distribution:", dist, "\n")
  cat("Type I error:", Type1_error[[i]], "\n")
  # cat("Statistics:\n")
  # print(test_statistic[[i]])
  cat("\n")
}


# Create density plots with theoretical t-density
plot_list <- vector("list", length(dist_sum))
for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  data <- data.frame(stat = test_statistic[[i]])
  
  # Calculate quartiles
  q2.5 <- qt(0.025, df = df)
  q97.5 <- qt(0.975, df = df)
  
  # Create plot
  p <- ggplot(data, aes(x = stat)) +
    geom_density(aes(color = "Empirical Density"), fill = "lightblue", alpha = 0.5) +
    geom_vline(aes(xintercept = q2.5, color = "2.5% Quartile"), linetype = "dashed") +
    geom_vline(aes(xintercept = q97.5, color = "97.5% Quartile"), linetype = "dashed") +
    stat_function(fun = dt, args = list(df = df), aes(color = "Theoretical t-Density"), linetype = "solid") +
    scale_color_manual(name = NULL,  # Remove legend title
                       values = c("Empirical Density" = "lightblue",
                                  "2.5% Quartile" = "red",
                                  "97.5% Quartile" = "blue",
                                  "Theoretical t-Density" = "green"),
                       labels = c("2.5% Quartile", "97.5% Quartile", "Empirical Density",  "Theoretical t-Density" )) +
    labs(title = paste(dist), x = "Test Statistic", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")  # Hide individual legends
  plot_list[[i]] <- p
}

# Create a dummy plot to extract the legend with proper colors and line types
dummy_data <- data.frame(x = c(1, 2, 3, 4), 
                         y = c(1, 2, 3, 4), 
                         color = factor(c("Empirical Density", "Theoretical t-Density", "2.5% Quartile", "97.5% Quartile")))

legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = color, linetype = color)) +
  geom_line(size = 1) +  # Use lines for legend
  scale_color_manual(values = c("Empirical Density" = "lightblue", 
                                "Theoretical t-Density" = "green",
                                "2.5% Quartile" = "red", 
                                "97.5% Quartile" = "blue"),
                     labels = c("2.5% Quartile", "97.5% Quartile", "Empirical Density", "Theoretical t-Density" )) +
  theme_void() +
  theme(legend.position = "right",  # Position the legend to the right
        legend.direction = "vertical", 
        legend.key = element_blank(),  # Remove box around legend labels
        legend.text = element_text(size = 10))  # Customize legend text

# Extract the legend from the dummy plot
legend <- get_legend(legend_plot)

# Combine plots into a 3x2 grid and place the legend on the right
combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = "v")
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(3, 0.4))  # Adjust relative widths to position the legend on the right

# Display the final plot
print(final_plot)


save.image(paste0("OneSample_Conditional_test_stat",".RData"))
