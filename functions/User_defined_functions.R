# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, dgof,  nortest, ggplot2, dplyr, tidyverse, LaplacesDemon, VGAM, patchwork, reshape2)

# ----------------------------------------------------
#   Generate data from different distribution 
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

# ---------------------------------------
# Define the functions to generate data 
# ---------------------------------------
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
generate_tests <- function(x, test){
  if(test == "KS"){
    output <- nortest::lillie.test(x)
  }
  else if(test == "SW"){
    output <- shapiro.test(x)
  }
  else if(test == "JB"){
    output <- tseries::jarque.bera.test(x)
  }
  else if(test == "DAP"){
    output <- moments::agostino.test(x)
  }
  else if(test == "AD"){
    output <- nortest::ad.test(x)
  }
  else if(test == "SF"){
    output <- nortest::sf.test(x)
  }
  else if(test == "CVM"){
    output <- nortest::cvm.test(x)
  }
  else if(test == "CHISQ"){
    output <- chisq.test(x)
  }else
  {
    stop("Unknown test name. Choose from: KS, SW, JB, DAP, AD, SF, CVM, CHISQ.")
  }
return(output)
}

# ------------------------------------
# Define simple calculation functions 
# -------------------------------------
# calculate one-sample  test statistic 
OneSample_test_statistic <- function(x) {
  return((mean(x)* sqrt(length(x))) / sd(x))
}
# calculate two-sample  test statistic 
TwoSample_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# compute the area under the curve 
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)) / (max(sample_size) - min(sample_size))
}

# ---------------------------------------------------------
# Define functions for performing various downstream tests 
# ---------------------------------------------------------

# -----------One Sample Permutation test ------
one_sample_permutation_test <- function(x, B){
  observe_stat <- OneSample_test_statistic(x)
  permuted_stat <- numeric(B)
  for (j in 1:B) {
    index <- sample(c(-1, 1), length(x), replace = TRUE)
    sample_data <- index * abs(x)
    permuted_stat[j] <- OneSample_test_statistic(sample_data)
  }
  
  return(mean(abs(permuted_stat) >= abs(observe_stat)))
}

# --------- Two Sample Permutation -----------
two_sample_permutation_test <- function(x, y, B) {
  observed <- TwoSample_test_statistic(x, y)
  combined <- c(x, y)
  perm_stats <- replicate(B, {
    permuted <- sample(combined)
    x_star <- permuted[1:length(x)]
    y_star <- permuted[(length(x)+1):(2*length(x))]
    TwoSample_test_statistic(x_star, y_star)
  })
  
  return(mean(abs(perm_stats) >= abs(observed)))
}

