# External functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Define the Generate_data function
Generate_data <- function(datagen.type, datagen.specify, n = NULL, dist = NULL, two_samples = FALSE, priors = NULL, x_weights = NULL, y_weights = NULL, ...) {
  
  # Case 1: Data generation via function
  # Check if sample size and distributions are specified
  if (datagen.type == 1) {
    if (is.null(n) || is.null(dist)) {
      stop("For datagen.type == 1, please specify both 'n' and 'dist' arguments.")
    }
    
    # Check if priors and weights are specified for a mixture distribution
    if (!is.null(priors) && length(priors) == 2 && !is.null(x_weights) && dist == "mixture") {
      if (length(x_weights) != 2) {
        stop("Please specify exactly two weights for x_weights for the mixture distribution.")
      }
      
      # Normalize weights for x to sum to 1
      x_weights <- x_weights / sum(x_weights)
      
      # Generate data for x using a mixture distribution with specified priors and weights
      x <- vector(length = n)
      for (i in 1:n) {
        component <- sample(1:2, size = 1, prob = x_weights)
        x[i] <- generate_data(1, priors[component])
      }
      #message("Sample x generated using specified mixture distribution.")
      
      # Generate data for y if two samples are specified
      if (two_samples) {
        if (is.null(y_weights)) {
          # If y_weights not provided, default to using x_weights
          y_weights <- x_weights
        } else {
          # Normalize weights for y to sum to 1
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
        #message("Sample y generated using specified mixture distribution.")
      } else {
        y <- NULL  # No second sample
      }
      
    } else {
      # Generate data for a single specified distribution if not using mixture
      x <- generate_data(n, dist)  # Generate data for x
      #message("Sample x generated using specified distribution.")
      
      if (two_samples) {
        y <- generate_data(n, dist)  # Generate data for y if two samples are specified
       # message("Sample y generated using specified distribution.")
      } else {
        y <- NULL  # No second sample
      }
    }
  }
  
  # Case 2: Data loading from a CSV file
  if (datagen.type == 2) {
    # Prompt user to select a CSV file
    file_path <- file.choose()  # Opens a dialog to select a file
    data <- read.csv(file_path, header = TRUE)  # Load the CSV data
    
    # Check if the file has enough columns based on the two_samples flag
    if (two_samples) {
      if (ncol(data) < 2) {
        stop("The CSV file should contain at least two columns for two samples (x and y).")
      }
      x <- data[[1]]  # First column for x
      y <- data[[2]]  # Second column for y
      message("Data loaded from CSV file with two samples.")
    } else {
      x <- data[[1]]  # Use the first column as data for x
      y <- NULL       # No second sample
      #message("Data loaded from CSV file with one sample.")
    }
  }
  
  # Return the data as a list for further processing
  return(list(x = x, y = y))
}

# generate data from a specified distribution
data1 <- Generate_data(datagen.type = 1, n = 10, dist = "Chi-Square", two_samples = FALSE)
data2 <- Generate_data(datagen.type = 1, n = 10, dist = "Chi-Square", two_samples = TRUE)
# Load two samples from a CSV:
data3 <- Generate_data(datagen.type = 2, two_samples = FALSE)
data4 <- Generate_data(datagen.type = 2, two_samples = TRUE)

# generate data from a mixture distribution
data_mixture <- Generate_data(
  datagen.type = 1, 
  n = n, 
  dist = "mixture", 
  priors = c("Exponential", "Standard Normal"), # specify priors for distributions
  x_weights = c(0, 10), # specify  weights for sample x
  y_weights = c(0, 10), # specify  weights for sample y
  two_samples = TRUE
)

# Example usage:
# To generate data from a mixture distribution with specified weights for x and y:
# Define Calculate_power function
Calculate_power <- function(d , alpha , N  , B = NULL, sample_size ) {
  powr_t <- numeric(length(sample_size))  # Initialize power results vector
  
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval_t <- numeric(N)  # Initialize vector to store p-values for each simulation
    
    for (i in 1:N) {
      # Generate data from a mixture distribution
      data_mixture <- Generate_data(
        datagen.type = 1, 
        n = n, 
        dist = "mixture", 
        priors = c("Exponential", "Standard Normal"),  # Specify priors for distributions
        x_weights = c(0, 10),                          # Weights for sample x
        y_weights = c(0, 10),                          # Weights for sample y
        two_samples = TRUE
      )
      
      # Perform downstream test using the Wilcox test
      pval_t[i] <- TwoSample_Power.test(data_mixture$x, data_mixture$y, "perm")
    }
    
    # Calculate power for current sample size
    powr_t[j] <- mean(pval_t < alpha)
  }
  
  return(powr_t)  # Return the vector of power values for each sample size
}

# Example usage of Calculate_power
power_results <- Calculate_power(d = 0.85, alpha = 0.05, N = 100, B = 100,  sample_size = c(10, 20, 30, 50))
power_results
