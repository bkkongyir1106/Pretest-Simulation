rm(list = ls())
# Install required packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest,gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Set working directory
setwd("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/Fall 2024/User Interface")
# Call in external functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# Define the Generate_data function
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
# Define Calculate_power function
Calculate_power <- function(alpha, N, twosamples = FALSE, dist, sample_size, test) {
  powr_t <- numeric(length(sample_size))
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval_t <- numeric(N)
    for (i in 1:N) {
      if (twosamples) {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = TRUE)
        pval_t[i] <- TwoSample_Power.test(data$x, data$y, test, alpha)
      } else {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = FALSE)
        pval_t[i] <- OneSample_Power.test(data$x, test, alpha)
      }
    }
    powr_t[j] <- mean(pval_t < alpha)
  }
  return(powr_t)
}
# Produce Power
set.seed(12345)
power_results <- Calculate_power(alpha = 0.05, N = 1000, twosamples = FALSE,
        dist = "Uniform", sample_size = c(10, 20, 30, 50, 100), test = "t_Wilcox")
power_results

#Calculate Type I Error
Calculate_TypeI_error <- function(alpha, N, twosamples = FALSE, dist, sample_size, test) {
  error <- numeric(length(sample_size))
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval_t <- numeric(N)
    for (i in 1:N) {
      if (twosamples) {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = TRUE)
        pval_t[i] <- TwoSample_TypeI_error_test(data$x, data$y, test, alpha)
      } else {
        data <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = FALSE)
        pval_t[i] <- OneSample_TypeI_test(data$x, test, alpha)
      }
    }
    error[j] <- mean(pval_t < alpha)
  }
  return(error)
}
# Produce Type I error
set.seed(12345)
TypeI.error_results <- Calculate_TypeI_error(alpha = 0.05, N = 1000, twosamples = TRUE,
                 dist = "Chi-Square", sample_size = c(10, 20, 30, 50, 100), test = "perm")
TypeI.error_results

# Bootstrap Method
data_load = Generate_data(datagen.type = 2, n = NULL, dist = NULL, two_samples = TRUE)
bootstrap_twosample.test(data_load$x, data_load$y)
