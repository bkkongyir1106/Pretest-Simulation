# Load dependencies and functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Type I error & power")

# Parallel setup
par_set <- function(cores_reserve = 2) {
  cores = parallel::detectCores()
  cores_use <- cores - cores_reserve
  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  return(cl)
}
close_cluster <- function(cl) parallel::stopCluster(cl)

# Bootstrap test for the mean
one_sample_bootstrap_test <- function(x, mu0 = 0, B = 1000) {
  T_obs <- mean(x) - mu0
  x_centered <- x - mean(x) + mu0
  boot_stats <- replicate(B, mean(sample(x_centered, replace = TRUE)) - mu0)
  p_val <- mean(abs(boot_stats) >= abs(T_obs))
  return(p_val)
}


# Parameters
{
  Nsim <- 1e5
  alpha <- 0.05
  B <- 1000
  effect_size <- 0.5  
  distributions <- c("Normal", "Uniform", "t", "Exponential", "Chi-Square", "LogNormal")
  sample_size <- c(8, 10, 15, 20, 25, 30, 50)
}

# Progress bar
{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size) * length(distributions)
  pb <- txtProgressBar(max = ntasks, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
}

# Simulation under H0 and H1
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM"),
                     .options.snow = opts) %:%
    foreach(dist = distributions) %dopar% {
      set.seed(12345)
      pval_t_H0 <- pval_wilcox_H0<-  pval_boot_H0 <- numeric(Nsim)
      pval_t_H1 <- pval_wilcox_H1 <- pval_boot_H1 <- numeric(Nsim)
      
      for (i in 1:Nsim) {
        
        x <- generate_data(n, dist)
        
        # Type I error
        pval_t_H0[i] <- t.test(x)$p.value
        pval_wilcox_H0[i] <- wilcox.test(x)$p.value
        pval_boot_H0[i] <- one_sample_bootstrap_test(x, mu0 = 0,  B)
        
        # power
        pval_t_H1[i] <- t.test(x + effect_size)$p.value
        pval_wilcox_H1[i] <- wilcox.test(x + effect_size)$p.value
        pval_boot_H1[i] <- one_sample_bootstrap_test(x + effect_size, mu0 = 0,  B)
        
      }
      
      list(
        # error
        error_t_test = mean(pval_t_H0 < alpha),
        error_wilcox_test = mean(pval_wilcox_H0 < alpha),
        error_bootstrap_test = mean(pval_boot_H0 < alpha),
        # power
        power_t_test = mean(pval_t_H1 < alpha),
        power_wilcox_test = mean(pval_wilcox_H1 < alpha),
        power_bootstrap_test = mean(pval_boot_H1 < alpha)
      )
    }
  
  close_cluster(my_cl)
  
  # store results
  
  errorvec <- numeric(length(sample_size) * length(distributions))
  TypeI_error_t.test <- TypeI_error_wilcox.test <- TypeI_error_bootstrap.test <- array(errorvec,
           dim = c(length(sample_size), length(distributions)),dimnames = list(sample_size, distributions))
  power_t.test <- power_wilcox.test <- power_bootstrap.test <- array(errorvec,
       dim = c(length(sample_size), length(distributions)),dimnames = list(sample_size, distributions))
  
  for (i in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      
      # Type I error
      TypeI_error_t.test[i, j] <- sim_out[[i]][[j]]$error_t_test
      TypeI_error_wilcox.test[i, j] <- sim_out[[i]][[j]]$error_wilcox_test
      TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
      
      # power
      power_t.test[i, j] <- sim_out[[i]][[j]]$power_t_test
      power_wilcox.test[i, j] <- sim_out[[i]][[j]]$power_wilcox_test
      power_bootstrap.test[i, j] <- sim_out[[i]][[j]]$power_bootstrap_test
    }
  }
})

# Print results
# Type I error
cat("Type I error Rates for t-test\n")
print(TypeI_error_t.test)

cat("Type I error Rates for Wilcoxon-test")
print(TypeI_error_wilcox.test)

cat("Type I error Rates for Bootstrap test\n")
print(TypeI_error_bootstrap.test)

# power
cat("power for t-test\n")
print(power_t.test)

cat("power for Wilcoxon-test")
print(power_wilcox.test)

cat("power for Bootstrap test\n")
print(power_bootstrap.test)

# Save results
save(Nsim, B, effect_size ,
     distributions, sample_size, 
     TypeI_error_t.test, 
     TypeI_error_wilcox.test,
     TypeI_error_bootstrap.test, 
     power_t.test, 
     power_wilcox.test,
     power_bootstrap.test, 
     file = "one_sample_test.RData"
)

# -----------------------------------------------
#                 plot Type I error
# -----------------------------------------------
pdf("one_sample_typeI_error.pdf", width = 10, height = 8)
# Layout: 2 rows x 3 columns for plots, 1 row for shared legend
layout_matrix <- matrix(1:7, nrow = 3, byrow = TRUE)
layout(mat = layout_matrix, heights = c(1, 1, 0.3))

# Plotting setup
par(mar = c(4, 4, 2, 1))
colors <- c("blue", "green", "red")
ltypes <- c(1, 2, 3)
pchs <- c(16, 17, 21)
tests <- c("t-test", "Wilcoxon test", "Bootstrap test")

for (j in seq_along(distributions)) {
  dist_name <- distributions[j]
  y_t <- TypeI_error_t.test[, j]
  y_w <- TypeI_error_wilcox.test[, j]
  y_b <- TypeI_error_bootstrap.test[, j]
  
  plot(sample_size, y_t, 
       type = "b", 
       col = colors[1], 
       lty = ltypes[1], 
       pch = pchs[1],
       lwd = 2,
       ylim = c(0, max(y_t, y_b, y_w, 0.2)),
       xlab = "Sample Size", 
       ylab = "Type I Error Rate",
       main = dist_name)
  
  lines(sample_size, y_w, type = "b", col = colors[2], lty = ltypes[1], pch = pchs[2], lwd = 2)
  lines(sample_size, y_b, type = "b", col = colors[3], lty = ltypes[1], pch = pchs[3], lwd = 2)
  abline(h = 0.05, col = "gray", lty = 3, lwd = 2)
}

# Legend (bottom row)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = tests,
       title = "Test Method", 
       col = colors, 
       lty = ltypes[1],
       lwd = 2,
       pch = pchs,
       horiz = TRUE, 
       bty = "n", 
       cex = 1.2)

dev.off()
# -----------------------------------------------
#                   Plot power
# -----------------------------------------------

pdf("one_sample_power.pdf", width = 10, height = 8)
# Layout: 2 rows x 3 columns for plots, 1 row for shared legend
layout_matrix <- matrix(1:7, nrow = 3, byrow = TRUE)
layout(mat = layout_matrix, heights = c(1, 1, 0.3))

# Plotting setup
par(mar = c(4, 4, 2, 1))
colors <- c("blue", "green", "red")
ltypes <- c(1, 2, 3)
pchs <- c(16, 17, 21)
tests <- c("t-test", "Wilcoxon test", "Bootstrap test")

for (j in seq_along(distributions)) {
  dist_name <- distributions[j]
  
  plot(sample_size, power_t.test[, j], 
       type = "b", 
       col = colors[1], 
       lty = ltypes[1], 
       pch = pchs[1],
       lwd = 2,
       ylim = c(0, 1),
       xlab = "Sample Size", 
       ylab = "Power",
       main = dist_name)
  
  lines(sample_size, power_wilcox.test[, j], type = "b", col = colors[2], lty = ltypes[1], pch = pchs[2], lwd = 2)
  lines(sample_size, power_bootstrap.test[, j], type = "b", col = colors[3], lty = ltypes[1], pch = pchs[3], lwd = 2)
  abline(h = 0.8, col = "gray", lty = 3)
}

# Legend (bottom row)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = tests,title = "Test Method", col = colors, lty = ltypes, pch = pchs,
       horiz = TRUE, bty = "n", cex = 1.2)

dev.off()
