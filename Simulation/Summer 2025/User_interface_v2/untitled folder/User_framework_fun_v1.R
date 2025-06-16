source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# --------------- loading user-defined functions ------------
# function to load user-defined function
load_user_test_function <- function(){
  file_path <- file.choose()
  env <- new.env()
  source(file = file_path, local = env)
  
  # get only function objects
  fun_names <- ls(env)
  funs_only <- fun_names[sapply(fun_names, function(x) is.function(env[[x]]))]
  
  # Error if no function is found
  if(length(funs_only) == 0){
    stop("No function found.")
  }
  
  # list all functions found
  cat("Functions found:\n")
  for (i in seq_along(funs_only)) {
    cat(i, ":", funs_only[i], "\n")
  }
  
  # prompt user to choose one
  if(length(funs_only) == 1){
    cat("Only one function found. Using it by default.\n")
    cat(funs_only, "function loaded successfully.\n")
    return(env[[funs_only[1]]])
  }else{
    repeat{
      choice <- as.integer(readline("Enter the number of the function to use: "))
      if(!is.na(choice) && choice >= 1 && choice <= length(funs_only)){
        cat(funs_only[choice], "function loaded successfully.\n")
        return(env[[funs_only[choice]]])
      }else{
        cat("Invalid input. Please enter a number between 1 and", length(funs_only), "\n")
      }
    }
  }
  
}

# ---------------- Extract and classify test results -------------
extract_results <- function(obj){
  out <- list()
  if(!is.null(obj$method)) out$method <- obj$method
  if(!is.null(obj$p.value)) out$p.value <- obj$p.value
  if(!is.null(obj$statistic)) out$statistic <- obj$statistic
  if(!is.null(obj$coefficients)) out$coefficients <- obj$coefficients
  
  # check if results contain a 'summary
  obj_classes <- class(obj)
  has_summary_method <- any(paste0("summary.", obj_classes) %in% methods("summary"))
  if(inherits(obj, "htest")){
    attr(out, "type") <- "test"
  }else if(has_summary_method){
    out$summary <- summary(obj)
    attr(out, "type") <- "model_summary"
  }else{
    out$result <- obj
    attr(out, "type") <- "raw"
  }
  return(out)
}

# ---------- Main wrapper - normality check & test run ----------
normality_check <- function(test_func, ..., data = NULL, test = "SW", alpha = 0.05) {
  args <- list(...)
  if (!"data" %in% names(args) && !is.null(data)) {
    args$data <- data
  }
  
  # Run the user's test function to prepare residuals if needed
  result <- do.call(test_func, args)
  
  # Prompt user for what to test
  cat("Which part of the test requires normality?\n")
  cat("Options: 'residuals', 'samples'\n")
  normality_target <- readline("Please enter your choice: ")
  
  if (!(normality_target %in% c("residuals", "samples"))) {
    stop("Invalid choice. Use either 'residuals' or 'samples'.")
  }
  
  # Run normality tests
  if (normality_target == "residuals") {
    normality_result <- generate_tests(residuals(result), test)
    cat("\n========Normality Test on Model's Residuals=====\n")
    cat("\n--------------------------------------\n")
    print(normality_result)
    if (is.null(normality_result$p.value)) stop("Residuals test did not return valid p-value.")
    norm_ok <- normality_result$p.value > alpha
  } else if (normality_target == "samples") {
    # Identify all numeric vectors in args for normality test
    sample_args <- args[sapply(args, is.numeric)]
    if (length(sample_args) == 0) stop("No numeric samples found in arguments.")
    
    normality_results <- lapply(sample_args, generate_tests, test = test)
    
    # Print all results
    for (i in seq_along(normality_results)) {
      cat("\n=====Normality test results for sample", i, "\n====")
      cat("\nSample", names(sample_args)[i], ":\n")
      cat("\n--------------------------------------\n")
      print(normality_results[[i]])
    }
    
    # Check if all p-values > alpha
    pvals <- sapply(normality_results, function(res) res$p.value)
    if (any(is.na(pvals))) stop("One or more normality tests failed to produce a p-value.")
    norm_ok <- all(pvals > alpha)
  }
  
  # Choose test based on normality
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
  
  # Extract and print result
  extracted <- extract_results(test_result)
  type <- attr(extracted, "type")
  
  if (type == "model_summary") {
    cat("\n========Model Summary=======\n")
    print(extracted$summary)
  } else if (type == "test") {
    cat("\n==========Test Results===========\n")
    if (!is.null(extracted$method)) cat("Method:", extracted$method, "\n")
    if (!is.null(extracted$statistic)) cat("Statistic:", extracted$statistic, "\n")
    if (!is.null(extracted$p.value)) cat("p-value:", extracted$p.value, "\n")
  } else {
    cat("\n=========Raw Results======\n")
    print(extracted$result)
  }
  
  invisible(list(
    normality_results = if (normality_target == "samples") normality_results else normality_result,
    final_test_result = test_result,
    extracted_summary = extracted
  ))
}

# ----------------------- Example Use ----------------------
# load functions
test_func <- load_user_test_function()
# Create dataset
set.seed(12345)
data <- list(a = generate_data(n = 50, dist = "Gamma"),
             b = generate_data(n = 50, dist = "Uniform"),
             c = generate_data(n = 50, dist = "Normal"),
             d = generate_data(n = 50, dist = "Normal"),
             x = generate_data(n = 50, dist = "Normal"),
             y = generate_data(n = 50, dist = "Exponential"),
             z = generate_data(n = 50, dist = "LogNormal"))

# Combine into a data frame for anova
group <- factor(rep(c("a","b" ,"x", "y", "z"), each = 50))
value <- c(data$a, data$b, data$x, data$y, data$z)
anova_data <- data.frame(group, value)

# ANOVA model
normality_check(test_func, formula = value ~ group, data = anova_data, test = "DAP")

# two sample test
normality_check(test_func, x = data$a, y = data$b,  test = "AD", alpha = 0.05)

# simple linear regression
normality_check(test_func, formula = data$z ~ data$a + data$b + data$x + data$y, data = data, test = "SW")
normality_check(test_func, formula = mpg ~ wt + hp + disp + drat, data = mtcars, test = "SW", alpha = 0.05)












