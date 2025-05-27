source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# --- step 1: Load user-defined functions from file ---
load_user_test_function <- function() {
  # browse for test function file
  file_path <- file.choose()
  # Load into an isolated environment.
  # This avoids overwriting existing variables/functions.
  env <- new.env()
  # Loads the R script into the new environment
  source(file_path, local = env)
  
  # Get only function objects
  fun_names <- ls(env)
  funs_only <- fun_names[sapply(fun_names, function(x) is.function(env[[x]]))]
  
  # Error if none found
  if (length(funs_only) == 0) {
    stop("No functions found.")
  }
  
  # List all functions found
  cat("Functions found:\n")
  for (i in seq_along(funs_only)) {
    cat(i, ":", funs_only[i], "\n")
  }
  
  # Prompt user to choose one
  if (length(funs_only) == 1) {
    cat("Only one function found. Using it by default.\n")
    cat("File loaded successfully.\n")
    return(env[[funs_only[1]]])
  } else {
    repeat {
      choice <- as.integer(readline("Enter the number of the function to use: "))
      if (!is.na(choice) && choice >= 1 && choice <= length(funs_only)) {
        cat("File loaded successfully.\n")
        return(env[[funs_only[choice]]])
      } else {
        cat("Invalid input. Please enter a number between 1 and", length(funs_only), "\n")
      }
    }
  }
}

# --- step 2: Extract and classify test results ---
extract_results <- function(obj) {
  out <- list()
  
  if (!is.null(obj$method)) out$method <- obj$method
  if (!is.null(obj$p.value)) out$p.value <- obj$p.value
  if (!is.null(obj$statistic)) out$statistic <- obj$statistic
  if (!is.null(obj$coefficients)) out$coefficients <- obj$coefficients
  
  # check if results object contains a 'summary'
  obj_classes <- class(obj)
  has_summary_method <- any(paste0("summary.", obj_classes) %in% methods("summary"))
  
  if (inherits(obj, "htest")) {
    attr(out, "type") <- "test"
  } else if (has_summary_method) {
    out$summary <- summary(obj)
    attr(out, "type") <- "model_summary"
  } else {
    out$result <- obj
    attr(out, "type") <- "raw"
  }
  
  return(out)
}


# --- step 3: Main wrapper - normality check & test run ---
normality_check <- function(test_func, ..., data = NULL, test = "SW", alpha = 0.05) {
  args <- list(...)
  if (!"data" %in% names(args) && !is.null(data)) {
    args$data <- data
  }
  
  # Run the user's test function
  result <- do.call(test_func, args)
  
  # Ask user which component to check for normality
  cat("Which part of the test requires normality?\n")
  cat("Options: 'residuals', 'x', 'y', 'diff', 'both'\n")
  normality_target <- readline("Enter your choice: ")
  
  valid_choices <- c("residuals", "x", "y", "diff", "both")
  if (!normality_target %in% valid_choices) {
    stop("Invalid choice. Use one of: ", paste(valid_choices, collapse = ", "))
  }
  
  # Run normality test
  normality_results <- switch(normality_target,
                              "residuals" = {
                                generate_tests(residuals(result), test)
                              },
                              "x" = generate_tests(args$x, test),
                              "y" = generate_tests(args$y, test),
                              "diff" = generate_tests(args$x - args$y, test),
                              "both" = list(
                                x_test = generate_tests(args$x, test),
                                y_test = generate_tests(args$y, test)
                              )
  )
  
  # Display and interpret normality test result
  cat("\n---", test, "Normality Test Result ---\n")
  if (normality_target == "both") {
    
    cat("Group x:\n"); print(normality_results$x_test)
    cat("Group y:\n"); print(normality_results$y_test)
    
    norm_ok <- normality_results$x_test$p.value > alpha &&
      normality_results$y_test$p.value > alpha
    
  } else {
    print(normality_results)
    if (!is.null(normality_results$p.value)) {
      
      norm_ok <- normality_results$p.value > alpha
      
    } else {
      stop("Normality test did not return a valid p-value.")
    }
  }
  
  # Perform downstream test
  if (norm_ok) {
    cat("\nNormality satisfied. Running parametric test...\n")
    test_result <- do.call(test_func, args)
  } else {
    cat("\nNormality violated.\n")
    cat("Please load a nonparametric alternative test function.\n")
    alt_func <- load_user_test_function()
    cat("Running nonparametric test...\n")
    test_result <- do.call(alt_func, args)
  }
  
  # Extract and classify test results
  extracted <- extract_results(test_result)
  
  # Display based on classification
  type <- attr(extracted, "type")
  if (type == "model_summary") {
    cat("\n--- Model Summary ---\n")
    print(extracted$summary)
  } else if (type == "test") {
    cat("\n--- Test Result ---\n")
    if (!is.null(extracted$method)) cat("Method:", extracted$method, "\n")
    if (!is.null(extracted$statistic)) cat("Statistic:", extracted$statistic, "\n")
    if (!is.null(extracted$p.value)) cat("p-value:", extracted$p.value, "\n")
  } else {
    cat("\n--- Raw Result ---\n")
    print(extracted$result)
  }
  
  invisible(list(
    normality_result = normality_results,
    final_test_result = test_result,
    extracted_summary = extracted
  ))
}

# load functions
test_func <- load_user_test_function()
# Run test
data <- list(x = generate_data(n = 10, dist = "Normal"),
             y = generate_data(n = 10, dist = "Normal"),
             z = generate_data(n = 10, dist = "Exponential"))

normality_check(test_func, x = data$x,  test = "AD", alpha = 0.05)

normality_check(test_func, formula = data$z ~ data$x + data$y, data = data, test = "SW")
normality_check(test_func, formula = mpg ~ wt + hp + disp + drat, data = mtcars, test = "SW", alpha = 0.05)
