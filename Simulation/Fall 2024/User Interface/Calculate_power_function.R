# Call in External functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Define the Generate_data function
Generate_data <- function(datagen.type, datagen.specify, n = NULL, dist = NULL, two_samples = FALSE, priors = NULL, x_weights = NULL, y_weights = NULL, ...) {
  
  # Case 1: Data generation via function
  if (datagen.type == 1) {
    if (is.null(n) || is.null(dist)) {
      stop("For datagen.type == 1, please specify both 'n' and 'dist' arguments.")
    }
    
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
      
      # Generate data for y if two samples are specified
      if (two_samples) {
        if (is.null(y_weights)) {
          # If y_weights not provided, default to using x_weights
          y_weights <- x_weights
        } else {
          if (length(y_weights) != 2) {
            stop("Please specify exactly two weights for y_weights for the mixture distribution.")
          }
          # Normalize weights for y to sum to 1
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
      # Generate data for a single specified distribution if not using mixture
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
  
  # Return the data as a list for further processing
  return(list(x = x, y = y))
}

# Calculate_power
Calculate_power <- function(d, alpha, N, B = NULL, sample_size, 
                            datagen.type = 1, dist = NULL, two_samples = FALSE, 
                            priors = NULL, x_weights = NULL, y_weights = NULL) {
  powr_t <- numeric(length(sample_size))  # Initialize power results vector
  
  for (j in seq_along(sample_size)) {
    n <- sample_size[j]
    pval_t <- numeric(N)  # Initialize vector to store p-values for each simulation
    
    for (i in 1:N) {
      # Generate data dynamically based on user input
      data_generated <- Generate_data(
        datagen.type = datagen.type, 
        n = n, 
        dist = dist, 
        priors = priors, 
        x_weights = x_weights, 
        y_weights = y_weights, 
        two_samples = two_samples
      )
      
      # Perform the appropriate test based on the number of samples generated
      if (is.null(data_generated$y)) {
        # One-sample test
        pval_t[i] <- OneSample_Power.test(data_generated$x, "t_Wilcox")
      } else {
        # Two-sample test
        pval_t[i] <- TwoSample_Power.test(data_generated$x, data_generated$y, "t_Wilcox")
      }
    }
    
    # Calculate power for current sample size
    powr_t[j] <- mean(pval_t < alpha)
  }
  
  return(powr_t)  # Return the vector of power values for each sample size
}

# Example usage of Calculate_power for different types of data generation:

# Case 1: Generate data for a single sample from a specified distribution
power_results_single_sample <- Calculate_power(
  d = 0.85, alpha = 0.05, N = 1000, 
  sample_size = c(10, 20, 30, 50, 100), 
  datagen.type = 1, dist = "Exponential", two_samples = FALSE
)

# Case 2: Generate data for two samples from a specified distribution
power_results_two_sample <- Calculate_power(
  d = 0.85, alpha = 0.05, N = 1000, 
  sample_size = c(10, 20, 30, 50, 100), 
  datagen.type = 1, dist = "Exponential", two_samples = TRUE
)

# Case 3: Generate data from a mixture distribution with specified priors and weights
power_results_mixture <- Calculate_power(
  d = 0.85, alpha = 0.05, N = 1000, 
  sample_size = c(10, 20, 30, 50, 100), 
  datagen.type = 1, dist = "mixture", two_samples = TRUE, 
  priors = c("Exponential", "t"), 
  x_weights = c(0.3, 0.7), y_weights = c(0.5, 0.5)
)

# We will skip this case for now.
# #Case 4: Load data from a CSV file (assuming datagen.type = 2)
# power_results_from_file <- Calculate_power(
#   d = 0.85, alpha = 0.05, N = 100,
#   sample_size = c(10, 20, 30, 50),
#   datagen.type = 2, two_samples = TRUE
# )

# Print results for each case
power_results_single_sample
power_results_two_sample
power_results_mixture
