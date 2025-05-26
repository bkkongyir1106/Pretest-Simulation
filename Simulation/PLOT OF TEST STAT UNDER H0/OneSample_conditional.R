rm(list = ls())

# Load libraries
library(ggplot2)
library(gridExtra)
library(cowplot)

# Set working directory and load custom functions
setwd("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/PLOT OF TEST STAT UNDER H0")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Parameters
dist_sum <- c("Normal", "Exponential",  "LogNormal")
set.seed(12345)
N <- 1e4
alpha <- 0.05
d <- 0.5
n <- 20
df <- n - 1

# Storage
test_statistic <- vector("list", length(dist_sum))
Type1_error <- vector("list", length(dist_sum))

# Main simulation loop
for (i in seq_along(dist_sum)) {
  dist <- dist_sum[i]
  Nsim <- 0
  Npassed <- 0
  stat <- numeric(N)
  
  while (Npassed < N) {
    x <- generate_data(n, dist)
    Nsim <- Nsim + 1
    if (shapiro.test(x)$p.value > alpha) {
      Npassed <- Npassed + 1
      stat[Npassed] <- t.test(x)$statistic
    }
  }
  
  test_statistic[[i]] <- stat
  left_norm <- mean(stat < qt(0.025, df = df))
  right_norm <- mean(stat > qt(0.975, df = df))
  Type1_error[[i]] <- left_norm + right_norm
}

# Display Type I errors
for (i in seq_along(dist_sum)) {
  cat("Distribution:", dist_sum[i], "\n")
  cat("Type I error:", Type1_error[[i]], "\n\n")
}

# Create density plots
plot_list <- vector("list", length(dist_sum))
for (i in seq_along(dist_sum)) {
  data <- data.frame(stat = test_statistic[[i]])
  dist <- dist_sum[i]
  
  p <- ggplot(data, aes(x = stat)) +
    geom_density(aes(color = "Empirical Density"), fill = "lightblue", alpha = 0.5) +
    geom_vline(aes(xintercept = qt(0.025, df = df), color = "2.5% Quartile"), linetype = "dashed") +
    geom_vline(aes(xintercept = qt(0.975, df = df), color = "97.5% Quartile"), linetype = "dashed") +
    stat_function(fun = dt, args = list(df = df), aes(color = "Theoretical t-Density")) +
    scale_color_manual(
      name = "Legend",
      values = c(
        "Empirical Density" = "lightblue",
        "2.5% Quartile" = "red",
        "97.5% Quartile" = "blue",
        "Theoretical t-Density" = "green"
      ),
      breaks = c("Empirical Density", "Theoretical t-Density", "2.5% Quartile", "97.5% Quartile")
    ) +
    labs(title = dist, x = "Test Statistic", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")
  
  plot_list[[i]] <- p
}

# Dummy plot to extract the legend
dummy_data <- data.frame(x = 1:4, y = 1:4, group = factor(c("Empirical Density", "Theoretical t-Density", "2.5% Quartile", "97.5% Quartile")))
legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1) +
  scale_color_manual(
    name = "Legend",
    values = c(
      "Empirical Density" = "lightblue",
      "Theoretical t-Density" = "green",
      "2.5% Quartile" = "red",
      "97.5% Quartile" = "blue"
    )
  ) +
  scale_linetype_manual(
    name = "Legend",
    values = c(
      "Empirical Density" = "solid",
      "Theoretical t-Density" = "solid",
      "2.5% Quartile" = "dashed",
      "97.5% Quartile" = "dashed"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  )

legend <- get_legend(legend_plot)

# Combine density plots in a 3x2 layout with the legend on the right
combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = "v")
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(3, 0.4))

# Show final plot
print(final_plot)

# Save session
save.image("OneSample_Conditional_test_stat.RData")


library(ggplot2)
library(cowplot)

# Create density plots function
generate_density_plots <- function(test_statistic, df) {
  plots <- vector("list", length(test_statistic))
  for (i in seq_along(test_statistic)) {
    stat_data <- data.frame(stat = test_statistic[[i]])
    q_low <- qt(0.025, df)
    q_high <- qt(0.975, df)
    
    plots[[i]] <- ggplot(stat_data, aes(x = stat)) +
      geom_density(aes(color = "Empirical Density"), fill = "lightblue", alpha = 0.5) +
      geom_vline(aes(xintercept = q_low, color = "2.5% Quantile"), linetype = "dashed") +
      geom_vline(aes(xintercept = q_high, color = "97.5% Quantile"), linetype = "dashed") +
      stat_function(fun = dt, args = list(df = df), aes(color = "Theoretical t-Density")) +
      scale_color_manual(
        name = "Legend",
        values = c(
          "Empirical Density" = "lightblue",
          "2.5% Quantile" = "red",
          "97.5% Quantile" = "blue",
          "Theoretical t-Density" = "green"
        )
      ) +
      labs(x = "Test Statistic", y = "Density") +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
      )
  }
  return(plots)
}

# Load all datasets and create plot sets
df <- 19
load("OneSample_test_stat.RData")
plots_os <- generate_density_plots(test_statistic, df)

load("OneSample_Conditional_test_stat.RData")
plots_os_cond <- generate_density_plots(test_statistic, df)

load("TwoSample_test_stat.RData")
plots_ts <- generate_density_plots(test_statistic, df)

load("TwoSample_Conditional_test_stat.RData")
plots_ts_cond <- generate_density_plots(test_statistic, df)

# Combine all plots into a list
all_plots <- c(plots_os, plots_os_cond, plots_ts, plots_ts_cond)

# Column and row labels
col_titles <- c("Normal", "Exponential",  "LogNormal")
row_labels <- c("One-Sample (Unconditional)", "One-Sample (Conditional)",
                "Two-Sample (Unconditional)", "Two-Sample (Conditional)")

# Create title grobs (column headers)
title_grobs <- lapply(col_titles, function(title) {
  ggdraw() + draw_label(title, fontface = "bold", size = 10, hjust = 0.5)
})
title_row <- plot_grid(plotlist = title_grobs, ncol = 6)

# Build each row with 6 plots + row label on right
make_row <- function(plots, label) {
  row <- plot_grid(plotlist = plots, ncol = 6)
  labeled_row <- plot_grid(row, ggdraw() + draw_label(label, size = 10, angle = 0, fontface = "bold"), 
                           ncol = 2, rel_widths = c(1, 0.07))
  return(labeled_row)
}

# Build all 4 rows
row1 <- make_row(plots_os, row_labels[1])
row2 <- make_row(plots_os_cond, row_labels[2])
row3 <- make_row(plots_ts, row_labels[3])
row4 <- make_row(plots_ts_cond, row_labels[4])

# Combine rows vertically
main_plot <- plot_grid(row1, row2, row3, row4, ncol = 1)

# Add column titles to top
full_plot_with_titles <- plot_grid(title_row, main_plot, ncol = 1, rel_heights = c(0.06, 1))

# Create and extract legend
dummy_data <- data.frame(
  x = 1:4,
  y = 1:4,
  group = factor(c("Empirical Density", "Theoretical t-Density", "2.5% Quantile", "97.5% Quantile"))
)
legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1) +
  scale_color_manual(
    name = "Legend",
    values = c(
      "Empirical Density" = "lightblue",
      "Theoretical t-Density" = "green",
      "2.5% Quantile" = "red",
      "97.5% Quantile" = "blue"
    )
  ) +
  scale_linetype_manual(
    name = "Legend",
    values = c(
      "Empirical Density" = "solid",
      "Theoretical t-Density" = "solid",
      "2.5% Quantile" = "dashed",
      "97.5% Quantile" = "dashed"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
legend <- get_legend(legend_plot)

# Combine everything into final layout
final_plot <- plot_grid(full_plot_with_titles, legend, ncol = 1, rel_heights = c(1, 0.08))

# Display the plot
print(final_plot)

