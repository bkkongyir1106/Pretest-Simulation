# Load necessary libraries
rm(list = ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/PLOT OF TEST STAT UNDER H0")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

dist_sum <- c("t", "Chi-Square", "LogNormal")
set.seed(12345)
N <- 1e4
alpha <- 0.05
n <- 10
df <- n - 1

# Initialize dataframe to store statistics from both simulations
test_stat_df <- data.frame()

# First Simulation (Same as before)
for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  stat <- numeric(N)  # Preallocate vector for statistics
  for (j in 1:N) {
    x <- generate_data(n, dist)
    stat[j] <- t.test(x)$statistic
  }
  
  # Combine the current distribution's statistics into the dataframe
  temp_df <- data.frame(Distribution = dist, Simulation = "Sim1", TestStatistic = stat)
  test_stat_df <- rbind(test_stat_df, temp_df)
}

# Second Simulation with Shapiro-Wilk pre-test
for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  Nsim2 <- 0
  Npassed2 <- 0
  stat2 <- numeric(N)  # Preallocate vector for statistics
  
  while (Npassed2 < N) {
    x <- generate_data(n, dist)
    Nsim2 <- Nsim2 + 1
    if (shapiro.test(x)$p.value > alpha) {  # Only accept if normality holds
      Npassed2 <- Npassed2 + 1
      stat2[Npassed2] <- t.test(x)$statistic
    }
  }
  
  # Combine the current distribution's statistics from the second simulation into the dataframe
  temp_df2 <- data.frame(Distribution = dist, Simulation = "Sim2", TestStatistic = stat2)
  test_stat_df <- rbind(test_stat_df, temp_df2)
}

# Create the density plots for each simulation
plot_list_sim1 <- vector("list", length(dist_sum))
plot_list_sim2 <- vector("list", length(dist_sum))

for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  
  # Filter data for the first simulation and the current distribution
  data_sim1 <- test_stat_df %>% filter(Simulation == "Sim1" & Distribution == dist)
  data_sim2 <- test_stat_df %>% filter(Simulation == "Sim2" & Distribution == dist)
  
  # Quartiles for t-distribution
  q2.5 <- qt(0.025, df = df)
  q97.5 <- qt(0.975, df = df)
  
  # Plot for Simulation 1
  p_sim1 <- ggplot(data_sim1, aes(x = TestStatistic)) +
    geom_density(aes(color = "Empirical Density"), fill = "lightblue", alpha = 0.5) +
    geom_vline(aes(xintercept = q2.5, color = "2.5% Quartile"), linetype = "dashed") +
    geom_vline(aes(xintercept = q97.5, color = "97.5% Quartile"), linetype = "dashed") +
    stat_function(fun = dt, args = list(df = df), aes(color = "Theoretical t-Density"), linetype = "solid") +
    scale_color_manual(name = "Legend",
                       values = c("Empirical Density" = "lightblue",
                                  "2.5% Quartile" = "green",
                                  "97.5% Quartile" = "blue",
                                  "Theoretical t-Density" = "red"),
                       labels = c("2.5% Quartile", "97.5% Quartile", "Empirical Density",  "Theoretical t-Density")) +
    labs(title = paste(dist, "- Unconditional"), x = "Test Statistic", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")  # Legend will be placed outside the plot
  
  plot_list_sim1[[i]] <- p_sim1
  
  # Plot for Simulation 2
  p_sim2 <- ggplot(data_sim2, aes(x = TestStatistic)) +
    geom_density(aes(color = "Empirical Density"), fill = "lightblue", alpha = 0.5) +
    geom_vline(aes(xintercept = q2.5, color = "2.5% Quartile"), linetype = "dashed") +
    geom_vline(aes(xintercept = q97.5, color = "97.5% Quartile"), linetype = "dashed") +
    stat_function(fun = dt, args = list(df = df), aes(color = "Theoretical t-Density"), linetype = "solid") +
    scale_color_manual(name = "Legend",
                       values = c("Empirical Density" = "lightblue",
                                  "2.5% Quartile" = "green",
                                  "97.5% Quartile" = "blue",
                                  "Theoretical t-Density" = "red"),
                       labels = c("2.5% Quartile", "97.5% Quartile", "Empirical Density",  "Theoretical t-Density")) +
    labs(title = paste(dist, "-Conditional"), x = "Test Statistic", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")  # Legend will be placed outside the plot
  
  plot_list_sim2[[i]] <- p_sim2
}

# Create a dummy plot to extract the legend with proper colors and line types (no boxes)
dummy_data <- data.frame(x = c(1, 2, 3, 4), 
                         y = c(1, 2, 3, 4), 
                         color = factor(c("Empirical Density", "Theoretical t-Density", "2.5% Quartile", "97.5% Quartile")))

legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = color, linetype = color)) +
  geom_line(size = 1) +  # Use lines for legend
  scale_color_manual(values = c("Empirical Density" = "lightblue", 
                                "Theoretical t-Density" = "red",
                                "2.5% Quartile" = "green", 
                                "97.5% Quartile" = "blue"),
                     labels = c("2.5% Quartile", "97.5% Quartile", "Empirical Density", "Theoretical t-Density")) +
  theme_void() +  # Empty plot background
  theme(legend.position = "right",  # Place legend on the right
        legend.direction = "vertical", 
        legend.key = element_blank(),  # Remove boxes around legend labels
        legend.text = element_text(size = 10))  # Customize legend text

# Extract the legend from the dummy plot
legend <- get_legend(legend_plot)

# Arrange plots in a 2x3 grid (Sim1 on top row, Sim2 on bottom row)
combined_plots <- plot_grid(
  plot_grid(plotlist = plot_list_sim1, ncol = 3, align = "v"),
  plot_grid(plotlist = plot_list_sim2, ncol = 3, align = "v"),
  ncol = 1
)

# Final plot with legend on the right
final_plot <- plot_grid(combined_plots, legend, ncol = 2, rel_widths = c(3, 0.4))

# Display the final plot
print(final_plot)

# Save the final plot if needed
# ggsave("density_plots_simulation_comparison_no_box_legend.png", final_plot, width = 16, height = 8)
