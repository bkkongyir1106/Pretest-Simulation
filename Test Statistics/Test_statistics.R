source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Set working directory
setwd("/Users/benedictkongyir/Library/Mobile Documents/com~apple~CloudDocs/PhD Thesis/Test Statistics")

# Simulation setup
set.seed(12345)
distribution <- c("Normal", "Exponential", "LogNormal")
N <- 1e5
alpha <- 0.05
n <- 10
df <- n - 1

# -------------------------
# One-sample t-test
# -------------------------
tstat_1s_uncond <- vector("list", length(distribution))
tstat_1s_cond   <- vector("list", length(distribution))

for (i in seq_along(distribution)) {
  dist <- distribution[i]
  
  # Unconditional
  stat <- numeric(N)
  for (j in 1:N) {
    x <- generate_data(n, dist)
    stat[j] <- t.test(x)$statistic
  }
  tstat_1s_uncond[[i]] <- stat
  
  # Conditional
  stat_c <- numeric(N)
  count <- 0
  while (count < N) {
    x <- generate_data(n, dist)
    if (shapiro.test(x)$p.value > alpha) {
      count <- count + 1
      stat_c[count] <- t.test(x)$statistic
    }
  }
  tstat_1s_cond[[i]] <- stat_c
}

# -------------------------
# Two-sample t-test
# -------------------------
tstat_2s_uncond <- vector("list", length(distribution))
tstat_2s_cond   <- vector("list", length(distribution))

for (i in seq_along(distribution)) {
  dist <- distribution[i]
  
  # Unconditional
  stat <- numeric(N)
  for (j in 1:N) {
    x <- generate_data(n, dist)
    y <- generate_data(n, dist)
    
    stat[j] <- t.test(x, y)$statistic
  }
  tstat_2s_uncond[[i]] <- stat
  
  # Conditional
  stat_c <- numeric(N)
  count <- 0
  while (count < N) {
    x <- generate_data(n, dist)
    y <- generate_data(n, dist)
    if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
      count <- count + 1
      stat_c[count] <- t.test(x, y)$statistic
    }
  }
  tstat_2s_cond[[i]] <- stat_c
}

# -------------------------
# Plotting: LogNormal Case
# -------------------------
pdf("tstat_density_lognormal.pdf", width = 12, height = 9)

i <- which(distribution == "LogNormal")
crit_vals <- qt(c(0.025, 0.975), df = df)
x_vals <- seq(-4, 4, length.out = 1000)
dens_theoretical <- dt(x_vals, df = df)

# Densities
dens_1s_uncond <- density(tstat_1s_uncond[[i]])
dens_1s_cond   <- density(tstat_1s_cond[[i]])
dens_2s_uncond <- density(tstat_2s_uncond[[i]])
dens_2s_cond   <- density(tstat_2s_cond[[i]])

# Plotting layout
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# One-sample plot
plot(x_vals, dens_theoretical, type = "l", col = "red", lwd = 2,
     ylim = c(0, max(dens_1s_uncond$y, dens_1s_cond$y, dens_theoretical)),
     xlab = "Test statistic", ylab = "Density", 
     main = "One-sample statistic (LogNormal)")
lines(dens_1s_uncond, col = "orange", lwd = 2)
lines(dens_1s_cond, col = "steelblue", lwd = 2, lty = 2)
abline(v = crit_vals, lty = 3)

# Two-sample plot
plot(x_vals, dens_theoretical, type = "l", col = "red", lwd = 2,
     ylim = c(0, max(dens_2s_uncond$y, dens_2s_cond$y, dens_theoretical)),
     xlab = "Test statistic", ylab = "Density", 
     main = "Two-sample statistic (LogNormal)")
lines(dens_2s_uncond, col = "orange", lwd = 2)
lines(dens_2s_cond, col = "steelblue", lwd = 2, lty = 2)
abline(v = crit_vals, lty = 3)

# Shared legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", legend = c("Theoretical t", "Unconditional", "Conditional", "Critical values"),
       col = c("red", "orange", "steelblue", "black"),
       lty = c(1, 1, 2, 3), lwd = c(2, 2, 2, 1),
       horiz = TRUE, bty = "n", cex = 1.2)

dev.off()

# -----------------------------------------
# PDF output
pdf("tstat_density_all_distributions.pdf", width = 12, height = 9)

# Plot layout: each distribution gets one row (2 columns: 1s and 2s)
n_dist <- length(distribution)
par(mfrow = c(n_dist, 2), mar = c(4, 4, 3, 1), oma = c(2, 0, 3, 0))

# Critical values
crit_vals <- qt(c(0.025, 0.975), df = df)
x_vals <- seq(-4, 4, length.out = 1000)
dens_theoretical <- dt(x_vals, df = df)

# Loop over distributions
for (i in seq_along(distribution)) {
  dist_name <- distribution[i]
  
  # One-sample
  d1_1s <- density(tstat_1s_uncond[[i]])
  d2_1s <- density(tstat_1s_cond[[i]])
  
  plot(x_vals, dens_theoretical, type = "l", col = "red", lwd = 2,
       ylim = c(0, max(d1_1s$y, d2_1s$y, dens_theoretical)),
       xlab = "Test statistic", ylab = "Density",
       main = paste("One-sample (", dist_name, ")", sep = ""))
  lines(d1_1s, col = "orange", lwd = 2)
  lines(d2_1s, col = "steelblue", lwd = 2, lty = 2)
  abline(v = crit_vals, lty = 3)
  
  # Two-sample
  d1_2s <- density(tstat_2s_uncond[[i]])
  d2_2s <- density(tstat_2s_cond[[i]])
  
  plot(x_vals, dens_theoretical, type = "l", col = "red", lwd = 2,
       ylim = c(0, max(d1_2s$y, d2_2s$y, dens_theoretical)),
       xlab = "Test statistic", ylab = "Density",
       main = paste("Two-sample (", dist_name, ")", sep = ""))
  lines(d1_2s, col = "orange", lwd = 2)
  lines(d2_2s, col = "steelblue", lwd = 2, lty = 2)
  abline(v = crit_vals, lty = 3)
}

# Add shared legend at bottom
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", legend = c("Theoretical t", "Unconditional", "Conditional", "Critical values"),
       col = c("red", "orange", "steelblue", "black"),
       lty = c(1, 1, 2, 3), lwd = c(2, 2, 2, 1),
       horiz = TRUE, bty = "n", cex = 1.2)

dev.off()
