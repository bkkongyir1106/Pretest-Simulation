# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, dgof,  nortest, ggplot2, dplyr, tidyverse, LaplacesDemon, VGAM)

# ----------------------------------------------------
#      Generate data from different distribution 
# ---------------------------------------------------
generate_data <- function(n, dist){
  if(dist == "Normal"){ 
    x <- rnorm(n, mean = 0, sd = 1)
  }
  if(dist == "Chi-Square"){
    x <- (rchisq(n, df = 3) - 3)/sqrt(6) 
  }
  if(dist == "Gamma"){
    x <- (rgamma(n, shape = 3, rate = 0.1) - 30)/sqrt(300)
  }
  if(dist == "Exponential"){
    x <- rexp(n, rate = 1) - 1 
  }
  if(dist == "t"){
    x <- (rt(n, df = 7))/sqrt(7/5) 
  }
  if(dist == "Uniform"){
    x <- (runif(n, min = 0, max = 1) - 0.5)*sqrt(12) 
  }
  if(dist == "Laplace"){
    x <- rlaplace(n, location = 0, scale = 4)/sqrt(32)
  }
  if(dist == "Weibull"){
    x <- (rweibull(n, shape = 1, scale = 2) - 2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2))) 
  }
  if(dist == "LogNormal"){
    x <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1)) 
  }
  
  if(dist == "Contaminated"){
    br <- rbinom(n , size = 1 , prob = 0.75)
    sd_br <- sqrt(1 + br * 24)
    x <- rnorm(n, mean = 0, sd = sd_br)/sqrt(0.25 * 1 + 0.75*25)
  }
  if(dist == "Pareto"){
    shape <- 3
    mean_pareto <- shape / (shape - 1)
    sd_pareto <- sqrt(shape / (((shape - 1)^2)*(shape - 2)))
    x <- (rpareto(n, shape = shape) - mean_pareto)/sd_pareto
  }
  return(x)
}

# -------------------------------------------------------------------------
# Define the functions to generate data from either a specific, 
# a mixture  distribution of load from a source file
# -------------------------------------------------------------------------
Generate_data <- function(datagen.type, n = NULL, dist = NULL, two_samples = FALSE, 
                          priors = NULL, x_weights = NULL, y_weights = NULL, ...) {
  
  # generation via function
  if (datagen.type == 1) {
    if (!is.null(priors) && length(priors) == 2 && !is.null(x_weights) && dist == "mixture") {
      x_weights <- x_weights / sum(x_weights)
      # Generate data for x 
      x <- vector(length = n)
      for (i in 1:n) {
        component <- sample(1:2, size = 1, prob = x_weights)
        x[i] <- generate_data(1, priors[component])
      }
      if (two_samples) {
        y_weights <- y_weights / sum(y_weights)
        # Generate data for y
        y <- vector(length = n)
        for (i in 1:n) {
          component <- sample(1:2, size = 1, prob = y_weights)
          y[i] <- generate_data(1, priors[component])
        }
      } else {
        y <- NULL  
      }
    } else {
      x <- generate_data(n, dist)  
      if (two_samples) {
        y <- generate_data(n, dist)  
      } else {
        y <- NULL  
      }
    }
  }
  # loading from a CSV file
  if (datagen.type == 2) {
    # Prompt user to select a CSV file
    file_path <- file.choose()  
    data <- read.csv(file_path, header = TRUE)  
    if (two_samples) {
      if (ncol(data) < 2) {
        stop("The CSV file should contain at least two columns for two samples (x and y).")
      }
      x <- data[[1]]  
      y <- data[[2]]  
    } else {
      x <- data[[1]]  
      y <- NULL       
    }
  }
  
  return(list(x = x, y = y))
}

# ---------------------------------------------------------------------------
# Define function to perform different Normality tests
# ---------------------------------------------------------------------------
# generate_tests <- function(x, test){
#   if(test == "KS"){
#     output <- lillie.test(x)
#   }
#   if(test == "SW"){
#     output <- shapiro.test(x)
#   }
#   if(test == "JB"){
#     output <- jarque.test(x)
#   }
#   if(test == "DAP"){
#     output <- agostino.test(x)
#   }
#   if(test == "AD"){
#     output <- ad.test(x)
#   }
#   if(test == "SF"){
#     output <- sf.test(x)
#   }
#   if(test == "CVM"){
#     output <- cvm.test(x)
#   }
#   if(test == "CHISQ"){
#     output <- chisq.test(x)
#   }
# return(output)
# }

# ----------------------- Modified tests ---------------------------
generate_tests <- function(data, test) {
  switch(test,
         KS = {
           # Note: KS test with estimated parameters is not recommended, consider Lilliefors test
           # result <- ks.test(data, "pnorm", mean(data), sd(data))
           # list(statistic = result$statistic, p.value = result$p.value)
           result <- lillie.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         SW = {
           result <- shapiro.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         JB = {
           result <- jarque.bera.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         DAP = {
           result <- agostino.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         AD = {
           # Corrected to use nortest's Anderson-Darling test
           result <- nortest::ad.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         SF = {
           result <- sf.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         CVM = {
           # Corrected to use nortest's Cramér-von Mises test
           result <- nortest::cvm.test(data)
           list(statistic = result$statistic, p.value = result$p.value)
         },
         stop("Invalid test specified")
  )
  return(result)
}

# -------------------------------------------------------------------
# ------------------------------------------------------------------
# Define simple calculation functions 
# ------------------------------------------------------------------
# calculate one-sample  test statistic 
OneSample_test_statistic <- function(x) {
  return((mean(x)  * sqrt(length(x))) / sd(x))
}
# calculate two-sample  test statistic 
TwoSample_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# compute the area under the curve 
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, - 1) + tail(y, - 1)) / 2)) / (max(sample_size) - min(sample_size))
}

# ---------------------------------------------------------------------------
# Define functions for performing various downstream tests 
# ---------------------------------------------------------------------------

# ------------------------------
# One Sample Case 
# -------------------------------
OneSample.test <- function(x, effect_size, test,  alpha, B = NULL){
  if(test == "t"){
    pval <- t.test(x, mu = effect_size)$p.value
  }
  if(test == "Wilcox"){
    pval <- wilcox.test(x, mu = effect_size)$p.value
  }
  if(test == "t_Wilcox"){
    if(shapiro.test(x)$p.value > alpha){
      pval <- t.test(x, mu = effect_size)$p.value
    }else{
      pval <- wilcox.test(x, mu = effect_size)$p.value
    }
  }# one sample permutation test only valid for symmetric distributions
  if(test == "perm"){
    observe_stat <- OneSample_test_statistic(x + effect_size)
    permuted_stat <- numeric(B)
    for (j in 1:B) {
      index <- sample(c(-1, 1), length(x), replace = TRUE)
      sample_data <- index * abs(x) 
      permuted_stat[j] <- OneSample_test_statistic(sample_data + effect_size)
    }
    pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }

  return(pval)
}

# ------------------------------
# Two Sample Case 
# -------------------------------
TwoSample.test <- function(x, y, effect_size, test, alpha,  B = NULL){
  if(test == "t"){
    pval <- t.test(x, y, mu = effect_size)$p.value
  }
  if(test == "Wilcox"){
    pval <- wilcox.test(x, y, mu = effect_size)$p.value
  }
  if(test == "t_Wilcox"){
    if(shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha){
      pval <- t.test(x, y, mu =  effect_size)$p.value
    }else{
      pval <- wilcox.test(x, y, mu = effect_size)$p.value
    }
  }# two sample permutation test valid for all distributions
  if(test == "perm"){
    observe_stat <- TwoSample_test_statistic(x, y + effect_size)
    data <- c(x, y)
    permuted_stat <- numeric(B)
    for (j in 1:B) {
      sample_data <- sample(data)
      sample_x <- sample_data[1:length(x)]
      sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
      permuted_stat[j] <- TwoSample_test_statistic(sample_x, sample_y + effect_size) 
    }
    pval <- mean(abs(permuted_stat) >= abs(observe_stat))
  }
 
  return(pval)
}

# ==============================================================================
# -----------------------Confidence Intervals Approach-------------------------
# bootstrap function for Type I error/power of one-sample location 
#                        test through Confidence intervals
# -----------------------------------------------------------------------------
{
  bootstrap_one_sample <- function(x, effect_size, alpha, n_bootstrap, sample_size) {
    results <- numeric(length(sample_size))
    
    for(j in seq_along(sample_size)) {
      n <- sample_size[j]
      
      # Null distribution
      null_stats <- replicate(n_bootstrap, {
        OneSample_test_statistic(sample(x, n, replace = TRUE), effect_size = 0)
      })
      crit <- quantile(null_stats, c(alpha/2, 1 - alpha/2))
      
      # Alternative distribution
      alt_stats <- replicate(n_bootstrap, {
        shifted_data <- sample(x, n, replace = TRUE)
        OneSample_test_statistic(shifted_data, effect_size)
      })
      
      results[j] <- mean(alt_stats < crit[1] | alt_stats > crit[2])
    }
    
    return(results)
  }
}

# -------------------------------------------------------------------------
# Example usage:
# ------------------------------------------------------------------------

{
  # sample_sizes <- c(10, 20, 30, 40, 50)  # Your sample sizes
  # #set.seed(12345)
  # power_results <- bootstrap_one_sample(rnorm(50), effect_size = 0.5, alpha = 0.05, n_bootstrap = 1e3, sample_size = sample_sizes)
  # 
  # # Create a data frame for plotting
  # plot_data <- data.frame(
  #   SampleSize = sample_sizes,
  #   Power = power_results
  # )
  # 
  # # Generate the ggplot
  # ggplot(plot_data, aes(x = SampleSize, y = Power)) +
  #   geom_line(color = "#4E79A7", linewidth = 1) +
  #   geom_point(color = "#E15759", size = 3) +
  #   scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  #   labs(x = "Sample Size",
  #        y = "Statistical Power",
  #        title = "Power Analysis by Sample Size",
  #        subtitle = "Bootstrap estimation of statistical power") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(face = "bold", hjust = 0.5),
  #     plot.subtitle = element_text(hjust = 0.5),
  #     panel.grid.minor = element_blank()
  #   )
  # 
}

# -----------------------------------------------------------------------
# -------------- One-Sample Bootstrap function for p-values -------------
# -----------------------------------------------------------------------
bootstrap_one_sample_test <- function(x, effect_size, alpha, n_bootstrap) {
  # Calculate observed test statistic
  T_obs <- OneSample_test_statistic(x, effect_size = 0)
  
  # Generate null distribution
  null_stats <- replicate(n_bootstrap, {
    
    sampled_x <- sample(x, length(x), replace = TRUE) 
    
    OneSample_test_statistic(sampled_x, effect_size)
  })
  
  # Calculate critical values (two-tailed)
  crit <- quantile(null_stats, c(alpha/2, 1 - alpha/2))
  
  # Return binary rejection decision for p-value
  as.integer(T_obs < crit[1] | T_obs > crit[2])
}

# ==============================================================================
# -----------------------Confidence Intervals Approach-------------------------
# bootstrap function for two-sample location test through Confidence intervals
# -----------------------------------------------------------------------------

{
  bootstrap_two_sample <- function(x1, x2, effect_size, alpha, n_bootstrap, sample_size) {
    results <- numeric(length(sample_size))
    #combined <- c(x1, x2, effect_size)
    
    for(j in seq_along(sample_size)) {
      n <- sample_size[j]
      
      # Null distribution
      null_stats <- replicate(n_bootstrap, {
        s1 <- sample(x1, n, replace = TRUE)
        s2 <- sample(x2, n, replace = TRUE)
        TwoSample_test_statistic(s1, s2, effect_size)
      })
      crit <- quantile(null_stats, c(alpha/2, 1 - alpha/2))
      
      # Alternative distribution
      alt_stats <- replicate(n_bootstrap, {
        s1 <- sample(x1, n, replace = TRUE)
        s2 <- sample(x2, n, replace = TRUE)  
        TwoSample_test_statistic(s1, s2, effect_size)
      })
      
      results[j] <- mean(alt_stats < crit[1] | alt_stats > crit[2])
    }
    
    return(results)
  }
}

# -------------------------------------------------------------------------
# Example usage:
# ------------------------------------------------------------------------
{
  # sample_sizes <- c(10, 20, 30, 40, 50)  # Your sample sizes
  # set.seed(12345)
  # power_results <- bootstrap_two_sample(rnorm(50), rnorm(50),  effect_size = 0.5, alpha = 0.05, n_bootstrap = 1e5, sample_size = sample_sizes)
  # # Create a data frame for plotting
  # plot_data <- data.frame(
  #   SampleSize = sample_sizes,
  #   Power = power_results
  # )
  # 
  # # Generate the ggplot
  # ggplot(plot_data, aes(x = SampleSize, y = Power)) +
  #   geom_line(color = "#4E79A7", linewidth = 1) +
  #   geom_point(color = "#E15759", size = 3) +
  #   scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  #   labs(x = "Sample Size",
  #        y = "Statistical Power",
  #        title = "Power Analysis by Sample Size",
  #        subtitle = "Bootstrap estimation of statistical power") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(face = "bold", hjust = 0.5),
  #     plot.subtitle = element_text(hjust = 0.5),
  #     panel.grid.minor = element_blank()
  #   )
}

# -----------------------------------------------------------------------
# -------------- two-sample bootstrap function for p-values -------------
# -----------------------------------------------------------------------
bootstrap_two_sample_test <- function(x1, x2, effect_size,  alpha, n_bootstrap) {
  # Calculate observed test statistic
  T_obs <- TwoSample_test_statistic(x1, x2, effect_size)
  
  # Create combined null distribution
  combined <- c(x1, x2 )
  n1 <- length(x1)
  n2 <- length(x2)
  
  # Generate null distribution
  null_stats <- replicate(n_bootstrap, {
    s1 <- sample(combined, n1, replace = TRUE)
    s2 <- sample(combined, n2, replace = TRUE) 
    TwoSample_test_statistic(s1, s2, effect_size)
  })
  
  # Calculate critical values (two-tailed)
  crit <- quantile(null_stats, c(alpha/2, 1 - alpha/2))
  
  # Return binary rejection decision
  as.integer(T_obs < crit[1] | T_obs > crit[2])
}


# --------------------Addressing Selection Effects ------------------
# -------------------- Data Splitting -------------------------------
# -------------------------------------------------------------------
{
  # set.seed(123)
  # dist_sum <- c("Normal", "t", "Exponential", "Chi-Square", "LogNormal")
  # N     <- 1e3
  # alpha <- 0.05
  # n     <- 20
  # iteration  <- c()
  # type_I_error <- numeric(length(dist_sum))
  # 
  # for(j in seq_along(dist_sum)){
  #   dist <- dist_sum[j]
  #   iter  <- 0
  #   pvals <- numeric(0)
  #   while (length(pvals) < N) {
  #     iter <- iter + 1
  #     x <- generate_data(n, dist)
  #   
  #     if (shapiro.test(x[1:floor(length(x) / 2)])$p.value > alpha) {
  #       pvals <- c(pvals, t.test(x[(floor(length(x) / 2) + 1):length(x)])$p.value)
  #     }
  #   }
  #   type_I_error[j] <- mean(pvals < alpha)
  #   iteration[j] <- iter
  # }
  
  # # Print Results
  # cat("Type I error rate:", type_I_error, "\n")
  # cat("Total iterations performed:", iteration, "\n")
  
}



