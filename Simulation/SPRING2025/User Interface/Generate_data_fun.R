
#source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# Define the Generate_data function
Generate_data <- function(datagen.type, n = NULL, dist = NULL, two_samples = FALSE, priors = NULL, x_weights = NULL, y_weights = NULL, ...) {
  
  # Case 1: Data generation via function
  if (datagen.type == 1) {
    if (is.null(n) || is.null(dist)) {
      stop("For datagen.type == 1, please specify both 'n' and 'dist' arguments.")
    }
    if (!is.null(priors) && length(priors) == 2 && !is.null(x_weights) && dist == "mixture") {
      if (length(x_weights) != 2) {
        stop("Please specify exactly two weights for x_weights for the mixture distribution.")
      }
      x_weights <- x_weights / sum(x_weights)
      
      # Generate data for x using a mixture distribution with specified priors and weights
      x <- vector(length = n)
      for (i in 1:n) {
        component <- sample(1:2, size = 1, prob = x_weights)
        x[i] <- generate_data(1, priors[component])
      }
      if (two_samples) {
        if (is.null(y_weights)) {
          y_weights <- x_weights
        } else {
          if (length(y_weights) != 2) {
            stop("Please specify exactly two weights for y_weights for the mixture distribution.")
          }
          y_weights <- y_weights / sum(y_weights)
        }
        
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
  
  # Case 2: Data loading from a CSV file
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
      message("Data loaded from CSV file with two samples.")
    } else {
      x <- data[[1]]  
      y <- NULL       
    }
  }
  
  return(list(x = x, y = y))
}

# generate data from a specified distribution
data1 <- Generate_data(datagen.type = 1, n = 10, dist = "Chi-Square", two_samples = FALSE)
data2 <- Generate_data(datagen.type = 1, n = 10, dist = "Chi-Square", two_samples = TRUE)
# Load two samples from a CSV:
data3 <- Generate_data(datagen.type = 2, two_samples = FALSE)
data4 <- Generate_data(datagen.type = 2, two_samples = TRUE)

# generate data from a mixture distribution
data_mixture <- Generate_data(datagen.type = 1, n = 20, dist = "mixture", 
                priors = c("Exponential", "Standard Normal"), x_weights = c(0, 10), 
                y_weights = c(0, 10), two_samples = TRUE)

# Define Calculate_power function
Calculate_power <- function(alpha , N, sample_size ) {
  powr_t <- numeric(length(sample_size))  
  
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval_t <- numeric(N)  
    
    for (i in 1:N) {
      pval_t[i] <- TwoSample_Power.test(data1$x, data_mixture$y, "t")
    }
    powr_t[j] <- mean(pval_t < alpha)
  }
  
  return(powr_t)  
}

# Example usage of Calculate_power
power_results <- Calculate_power(alpha = 0.05, N = 100, sample_size = c(10, 20, 30, 50))
power_results
