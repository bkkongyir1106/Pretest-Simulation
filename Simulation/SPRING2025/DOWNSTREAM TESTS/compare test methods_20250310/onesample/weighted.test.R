# Load user-defined functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/onesample")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Define parallel cluster setup
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
  cat("Using", foreach::getDoParWorkers(), "cores\n")
  return(cl)
}

close_cluster <- function(cl) {
  parallel::stopCluster(cl)
}

# Simulation parameters
Nsim <- 1e3
sample_size <- c(8, 10, 15, 20, 25, 30)
distribution <- c("Normal", "Exponential", "LogNormal")
alpha <- 0.05
n_boot <- 1e3
effect_size <- 0.0  # Type I error case (null hypothesis is true)

# Setup parallelization
my_cl <- par_set(cores_reserve = 2)
ntasks <- length(sample_size) * length(distribution)
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run simulation
system.time({
  sim_out <- foreach(i = 1:ntasks,
                     .packages = c("LaplacesDemon", "VGAM", "snow"),
                     .options.snow = opts) %dopar% {
                       n <- sample_size[((i - 1) %% length(sample_size)) + 1]
                       dist <- distribution[((i - 1) %/% length(sample_size)) + 1]
                       set.seed(1234)
                       
                       # Initialize storage
                       pval.t_test <- numeric(Nsim)
                       pval.boot.test <- numeric(Nsim)
                       pval.adaptive.t.test <- c()
                       pval.adaptive.boot.test <- c()
                       pval.adaptive_split.t.test <- c()
                       pval.adaptive_split.boot.test <- c()
                       count_adaptive.t.test <- 0
                       count_adaptive.boot.test <- 0
                       count_adaptive_split.t.test <- 0
                       count_adaptive_split.boot.test <- 0
                       
                       for (j in 1:Nsim) {
                         x <- generate_data(n, dist)
                         
                         # Perform tests
                         pval.t_test[j] <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size)
                         pval.boot.test[j] <- one_sample_bootstrap_p_value(x, effect_size = effect_size, n_bootstrap = n_boot)
                         
                         # Adaptive test
                         if (shapiro.test(x)$p.value > alpha) {
                           count_adaptive.t.test <- count_adaptive.t.test + 1
                           pval.adaptive.t.test[count_adaptive.t.test] <- OneSample.test(x, test = "t", alpha = alpha, effect_size = effect_size)
                         } else {
                           count_adaptive.boot.test <- count_adaptive.boot.test + 1
                           pval.adaptive.boot.test[count_adaptive.boot.test] <- one_sample_bootstrap_p_value(x, effect_size = effect_size, n_bootstrap = n_boot)
                         }
                         
                         # Adaptive split test
                         split_point <- floor(n / 2)
                         first_half.x <- x[1:split_point]
                         second_half.x <- x[(split_point + 1):n]
                         
                         if (shapiro.test(first_half.x)$p.value > alpha) {
                           count_adaptive_split.t.test <- count_adaptive_split.t.test + 1
                           pval.adaptive_split.t.test[count_adaptive_split.t.test] <- OneSample.test(second_half.x, test = "t", alpha = alpha, effect_size = effect_size)
                         } else {
                           count_adaptive_split.boot.test <- count_adaptive_split.boot.test + 1
                           pval.adaptive_split.boot.test[count_adaptive_split.boot.test] <- one_sample_bootstrap_p_value(second_half.x, effect_size = effect_size, n_bootstrap = n_boot)
                         }
                       }
                       
                       # Aggregate Type I error results
                       Results <- list(
                         error.t.test = mean(pval.t_test < alpha),
                         error.boot.test = mean(pval.boot.test < alpha),
                         error.adaptive_boot.test = ((count_adaptive.t.test)/Nsim) * mean(pval.adaptive.t.test < alpha) +
                           ((count_adaptive.boot.test)/Nsim) * mean(pval.adaptive.boot.test < alpha),
                         error.adaptive_split.test = ((count_adaptive_split.t.test)/Nsim) * mean(pval.adaptive_split.t.test < alpha) +
                           ((count_adaptive_split.boot.test)/Nsim) * mean(pval.adaptive_split.boot.test < alpha)
                       )
                       return(Results)
                     }
})

close_cluster(my_cl)

# Initialize result matrices
powervec <- numeric(length(sample_size) * length(distribution))
TypeIerror.t.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.boot.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.adaptive_boot.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))
TypeIerror.adaptive_split.test <- array(powervec, dim = c(length(sample_size), length(distribution)), dimnames = list(sample_size, distribution))

# Populate final results
for (i in seq_along(sample_size)) {
  for (j in seq_along(distribution)) {
    index <- (j - 1) * length(sample_size) + i
    TypeIerror.t.test[i, j] <- sim_out[[index]]$error.t.test
    TypeIerror.boot.test[i, j] <- sim_out[[index]]$error.boot.test
    TypeIerror.adaptive_boot.test[i, j] <- sim_out[[index]]$error.adaptive_boot.test
    TypeIerror.adaptive_split.test[i, j] <- sim_out[[index]]$error.adaptive_split.test
  }
}

# Optional: Print results
print(TypeIerror.t.test)
print(TypeIerror.boot.test)
print(TypeIerror.adaptive_boot.test)
print(TypeIerror.adaptive_split.test)
