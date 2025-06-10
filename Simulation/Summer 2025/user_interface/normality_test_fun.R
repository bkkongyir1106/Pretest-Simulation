normality_test_for_residuals <- function(model, test = "SW"){
  normality_test <- generate_tests(residuals(model), test = test)
  cat("\n", test, "test for normality\n")
  return(list(
    test_statistic = normality_test$statistic,
    p_value        = normality_test$p.value
  ))

}

normality_test_for_samples <- function(samples, test = "SW") {
  results <- vector("list", length(samples))
  
  for (i in seq_along(samples)) {
    cat("\nSample", i, "-", test, "test for normality\n")
    results[[i]] <- generate_tests(samples[[i]], test = test)
  }
  
  return(results)
}

#Example
normality_test_for_samples(samples = list(rnorm(10), rexp(10), rchisq(10, df = 3)), test = "SW")




