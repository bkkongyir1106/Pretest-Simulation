#rm(list = ls())
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("/Users/benedictkongyir/Library/Mobile Documents/com~apple~CloudDocs/PhD Thesis/Normality_test_Methods")
# set up cores for parallel processing
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores()
  cores_use <- cores - cores_reserve
  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  foreach::getDoParWorkers()
  return(cl)
}

close_cluster <- function(cl) {
  parallel::stopCluster(cl) 
}
# ------------ set up simulation ----------------------
{
N <- 1e5
alpha <- 0.05
distributions <- c("Normal", "Uniform", "t", "Laplace", "Exponential", 
                "Chi-Square", "Gamma", "Weibull", "LogNormal", "Contaminated", "Pareto")
testvec <- c("KS", "SW", "JB", "DAP", "AD", "SF", "CVM")
sample_size <- c(8, 10, 15, 20, 25, 30, 50)
}

{
  my_cl <- par_set(cores_reserve = 2)
  ntasks <- length(sample_size) * length(testvec) * length(distributions)
  pb <- txtProgressBar(max = ntasks, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
}
# ----------------------------------------------------
## Perform simulation
system.time({
  sim_out <- foreach(n = sample_size,
                     .packages = c("LaplacesDemon", "VGAM", "moments", "nortest"),
                     .options.snow = opts) %:% 
              foreach(dist = distributions) %:% 
              foreach(test = testvec) %dopar% {   
      set.seed(12345)  
      pval <- numeric(N)
      for (i in 1:N) {
        x <- generate_data(n, dist)  
        pval[i] <- generate_tests(x, test)$p.value
      }
      powr <- mean(pval < alpha)
      list(powr = powr)
    }
  
  close_cluster(my_cl)        
  
  ## Output
  powervec <- numeric(length(sample_size) * length(distributions) * length(testvec))
  power_norm.test <- array(powervec, 
                           dim = c(length(sample_size), length(distributions), length(testvec)), 
                           dimnames = list(sample_size, distributions, testvec))
  for (t in seq_along(sample_size)) {
    for (j in seq_along(distributions)) {
      for (k in seq_along(testvec)) {
        power_norm.test[t, j, k] <- sim_out[[t]][[j]][[k]]$powr
      }
    }
  }
})

power_norm.test

save(
  N, 
  distributions,
  testvec,
  sample_size,
  power_norm.test,
  file = "Power_Normality_test.RData"
)


# Load simulation results
load("Power_Normality_test.RData")

pdf("Power_of_normality_test_methods.pdf", width = 10, height = 6)
# Select distributions
selected_distributions <- c("Uniform", "t", "Laplace", "Exponential", "Chi-Square", "LogNormal")
selected_indices <- match(selected_distributions, distributions)

# Define layout: 2 rows for plots + 1 short row for legend
layout(matrix(c(1:6, 7, 7, 7), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.2))

# Plotting symbols and colors
shapes <- 1:length(testvec)
colors <- rainbow(length(testvec))

# Set common margins
par(mar = c(4, 4, 2, 1))

# Plot each selected distribution
for (j in selected_indices) {
  plot(NULL, type = "n",
       xlim = range(sample_size),
       ylim = c(0, 1),
       xlab = "Sample Size",
       ylab = "Power",
       main = paste("Power for", distributions[j]))
  
  for (k in seq_along(testvec)) {
    yvals <- power_norm.test[, j, k]
    lines(sample_size, yvals, type = "b", pch = shapes[k], col = colors[k], lwd = 2)
  }
}

# Plot shared legend at the bottom
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = testvec,
       col = colors,
       pch = shapes,
       lty = 1,
       lwd = 2,
       horiz = TRUE,
       bty = "n",
       cex = 0.9,
       title = "Test")

dev.off()