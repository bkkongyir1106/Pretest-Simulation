source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# ------- Extract the part to test for normality ------
calculate_residuals <- function(model){
  return(resid(model))
}

calculate_samples <- function(samples) {
  if (is.data.frame(samples) || is.matrix(samples)) {
    return(as.list(as.data.frame(samples)))
  } else if (is.list(samples)) {
    return(samples)
  } else {
    return(list(samples))  
  }
}

#------------- perform normality test ---------------
normality_test <- function(data, test = "SW", alpha = 0.05) {
  pvals <- NULL
  
  # Case 1: Single numeric vector
  if (is.numeric(data) && is.atomic(data) && is.null(dim(data))) {
    pval <- generate_tests(data, test = test)$p.value
    cat("\np-value for", test, "test for normality\n")
    return(list(
      p_values = pval,
      normality_satisfied = pval > alpha
    ))
  }
  
  # Case 2: List of numeric vectors
  else if (is.list(data) && !is.data.frame(data)) {
    pvals <- sapply(data, function(sample) {
      if (!is.numeric(sample)) 
        stop("All elements in the list must be numeric vectors.")
      generate_tests(sample, test = test)$p.value
    })
    names(pvals) <- names(data) %||% paste0("Sample", seq_along(pvals))
  }
  
  # Case 3: Wide-format data frame (each column is a numeric sample)
  else if (is.data.frame(data) && all(sapply(data, is.numeric))) {
    sample_list <- as.list(data)
    pvals <- sapply(sample_list, function(x) generate_tests(x, test = test)$p.value)
    names(pvals) <- names(data)
  }
  
  # Case 4: Long-format data frame/matrix with 2 columns
  else if ((is.data.frame(data) || is.matrix(data)) && ncol(data) >= 2) {
    grouped_samples <- split(data[[1]], data[[2]])
    pvals <- sapply(grouped_samples, function(sample) {
      generate_tests(sample, test = test)$p.value
    })
  }
  
  # Otherwise unsupported
  else {
    stop("Unsupported input type: must be numeric vector, list of vectors, or 2-column data.frame/matrix.")
  }
  
  # Final result
  cat("\np-values for", test, "test for normality\n")
  return(list(
    p_values = pvals,
    normality_satisfied = all(pvals > alpha)
  ))
}

# Example
samples <- list(norm = rnorm(5), exp = rexp(5), chisq = rchisq(5, df = 3), t = rt(5, df = 3))
normality_test(samples, test = "SW", alpha = 0.05)

# linear regression
model <- lm(mpg ~ disp + hp + drat, data = mtcars)
samples <- residuals(model)

normality_test(samples, test = "AD", alpha = 0.05)

# anova
anova_model <- aov(mpg ~ as.factor(cyl), data = mtcars)
resid_aov <- residuals(anova_model)
normality_test(resid_aov, test = "SW", alpha = 0.05)
