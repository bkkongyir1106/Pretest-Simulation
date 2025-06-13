# Identify what needs to be tested for normality
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
# -------------------------------------------

# -------------------test for normality -------------

normality_test_from_function <- function(f, input, test = "SW") {
  
  result <- f(input)
  
  if (is.numeric(result) && is.atomic(result)) {
    return(generate_tests(result, test = test)$p.value)
  }
  
  if (is.list(result)) {
    pvals <- sapply(result, function(x) 
      generate_tests(x, test = test)$p.value)
    names(pvals) <- paste0("Sample", seq_along(pvals))
    return(pvals)
  }
  
  stop("Unsupported result type returned from function.")
}
# -------------------------------------------------
generate_data_function <- function(n, 
                                   para = list(beta0 = 0, beta1 = 1), 
                                   error_generator = rnorm, 
                                   error_args = list(mean = 0, sd = 6), 
                                   f = calculate_residuals) {
  # 1. Generate predictor
  x <- runif(n)
  z <- rexp(n)
  a <- rnorm(n)
  
  # 2. Generate error 
  error <- do.call(error_generator, c(list(n = n), error_args))
  
  # 3. Compute response
  y <- para$beta0 + para$beta1 * x + error
  
  # 4. Build data
  df <- data.frame(x = x, y = y, z = z, a = a)
  
  # 5. Fit model
  model <- lm(y ~ x + z, data = df)
  
  # 6.return part specified by f()
  return(f(df))
}
# --------------- Example ------------------
# return the specified components for normality
set.seed(12345)
generate_data_function(
  n = 10,
  para = list(beta0 = 3, beta1 = 7),
  error_generator = rnorm,
  error_args = list(mean = 50, sd = 15),
  f = calculate_samples
)
# -------- test for normality ----------
# Example 2: exponential errors,  residuals
normality_test_from_function(
  f = function(df) generate_data_function(
    n = 30,
    para = list(beta0 = 3, beta1 = 7),
    error_generator = rexp,
    error_args = list(rate = 3),
    f = calculate_samples
  ),
  input = NULL,
  test = "SW"
)
# ------------------------------------------------

# -----------------Downstream test ---------------
# two-sample t-test
two_ind_t_test <- function(x, y = NULL, mu = 0) {
  return(t.test(x = x, y = y, mu = mu, paired = FALSE))
}

two_dep_t_test <- function(x, y = NULL, mu = 0) {
  return(t.test(x = x, y = y, mu = mu, paired = TRUE))
}

one_sample_t_test <- function(x, mu = 0) {
  return(t.test(x = x, mu = mu))
}

# linear regression
simple_linear_regression <- function(formula, data) {
  return(summary(lm(formula, data = data)))
}

# Mann Whitney U test
Mann_Whitney_U_test <- function(){
  return(wilcox.test(x, y))
}

anova_test <- function(formula, data){
  return(aov(formula, data))
}
