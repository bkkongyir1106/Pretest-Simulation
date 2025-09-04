gen_data <- function(n = 20, effect_size = 0.0, sd = 1, dist = "Exponential", par = NULL
) {
  dist <- tolower(dist)
  x <- effect_size + sd * generate_data(n, dist)
  
  return(x)
}

# get parameters
get_parameters <- function(n, effect_size = 0.5, sd = 1, dist = "exponential", par = NULL, ...) {
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par
  )
}

# get normality object
fn_to_get_norm_obj <- function(data) {
  return(data)
}

# ds test function 1
fn_for_ds_test_1 <- function(data) {
  test_result <- t.test(data, mu = 0)
  return(list(p.value = test_result$p.value))
}

# ds test function 2: Sign test for the median
fn_for_ds_test_2 <- function(data, mu0 = 0) {
  x0 <- data - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  return(list(p.value = binom.test(signs, n_valid, p = 0.5)$p.value))
}
