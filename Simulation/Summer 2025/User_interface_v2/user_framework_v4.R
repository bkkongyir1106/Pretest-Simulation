# Load necessary utilities
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

# ---------------------------------------------------------------
# Load a user-defined function from a file
load_user_test_function <- function() {
  file_path <- file.choose()
  env <- new.env()
  source(file = file_path, local = env)
  
  fun_names <- ls(env)
  funs_only <- fun_names[sapply(fun_names, function(x) is.function(env[[x]]))]
  
  if (length(funs_only) == 0) {
    stop("No function found in the file.")
  }
  
  cat("Functions found:\n")
  for (i in seq_along(funs_only)) {
    cat(i, ":", funs_only[i], "\n")
  }
  
  if (length(funs_only) == 1) {
    cat("Only one function found. Using it by default.\n")
    return(env[[funs_only[1]]])
  } else {
    repeat {
      choice <- as.integer(readline("Enter the number of the function to use: "))
      if (!is.na(choice) && choice >= 1 && choice <= length(funs_only)) {
        cat(funs_only[choice], "function loaded successfully.\n")
        return(env[[funs_only[choice]]])
      } else {
        cat("Invalid input. Please enter a number between 1 and", length(funs_only), "\n")
      }
    }
  }
}

# ----------------------------------------------------------------
# Main function to perform normality test and subsequent analysis
normality_test <- function(test = "SW", alpha = 0.05) {
  # Step 1: Load data-generation function
  cat("Please load your data generation function")
  generate_function <- load_user_test_function()
  
  # Generate data
  generated_data <- generate_function()
  
  # Ask if normality should be tested on raw data
  response <- tolower(readline(prompt = "\nIs normality to be assessed on the raw data? (yes/no):\n"))
  
  if (response == "yes") {
    data <- generated_data
  } else {
    # Step 2: Load function that structures/extracts data for testing
    cat("\nChoose a function to extract part of model/data to test for normality\n")
    extractor_function <- load_user_test_function()
    
    # Extract target structure
    data <- extractor_function(generated_data)
  }
  
  # Initialize variables for results
  pvals <- NULL
  all_normality_satisfied <- NULL
  
  # Case 1: Single numeric vector
  if (is.numeric(data) && is.atomic(data) && is.null(dim(data))) {
    pval <- generate_tests(data, test = test)$p.value
    pvals <- pval
    all_normality_satisfied <- pval > alpha
  }
  
  # Case 2: List of numeric vectors
  else if (is.list(data) && !is.data.frame(data)) {
    pvals <- sapply(data, function(sample) {
      if (!is.numeric(sample))
        stop("All elements in the list must be numeric vectors.")
      generate_tests(sample, test = test)$p.value
    })
    names(pvals) <- names(data) %||% paste0("Sample", seq_along(pvals))
    all_normality_satisfied <- all(pvals > alpha)
  }
  
  # Case 3: Wide-format data frame
  else if (is.data.frame(data) && all(sapply(data, is.numeric))) {
    pvals <- sapply(as.list(data), function(x) generate_tests(x, test = test)$p.value)
    names(pvals) <- names(data)
    all_normality_satisfied <- all(pvals > alpha)
  }
  
  # Case 4: Long-format with group labels
  else if ((is.data.frame(data) || is.matrix(data)) && ncol(data) >= 2) {
    # We assume first column is group, second is value
    grouped_samples <- split(data[[2]], data[[1]])
    pvals <- sapply(grouped_samples, function(sample) {
      generate_tests(sample, test = test)$p.value
    })
    all_normality_satisfied <- all(pvals > alpha)
  }
  
  # Unsupported format
  else {
    stop("Unsupported input type: must be numeric vector, list of vectors, or grouped data.")
  }
  
  # Print normality results
  cat("\nNormality test results:\n")
  print(pvals)
  cat("Normality satisfied:", all_normality_satisfied, "\n")
  
  # Ask if user wants to run main statistical test
  run_main <- tolower(readline(prompt = "\nDo you want to run the main statistical test? (yes/no): \n"))
  
  if (run_main != "yes") {
    cat("\nExiting without running main test.\n")
    return(list(
      normality_results = list(p_values = pvals, 
                               normality_satisfied = all_normality_satisfied)
    ))
  }
  
  # Load main test function
  cat("\nChoose the function for the main statistical test\n")
  main_test_function <- load_user_test_function()
  
  # Run main test on original generated data
  cat("\nRunning main statistical test...\n")
  test_result <- main_test_function(generated_data)
  
  # Print the test results
  cat("\nMain Test Results:\n")
  print(test_result)
  
  # Return both sets of results
  return(list(
    normality_results = list(p_values = pvals, 
                             normality_satisfied = all_normality_satisfied),
    main_test_results = test_result
  ))
}

# Example use
set.seed(12345)
results <- normality_test(test = "AD", alpha = 0.05)
print(results)
