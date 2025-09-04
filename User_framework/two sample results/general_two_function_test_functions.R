# Data generation function
gen_data <- function(
    n1 = n,         
    n2 = n,         
    mean1 = 0.0,       
    effect_size = 0.0,     
    sd1 = 1,         
    sd2 = 1,         
    dist = "chi_square",
    par = NULL
) {
  
  group1 <- mean1 + sd1 * generate_data(n1, dist, par = par)
  group2 <- effect_size + sd2 * generate_data(n2, dist, par = par)
  
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
}

# parameters function
get_parameters <- function(n, ...) {
  defaults <- list(
    n1 = n,
    n2 = n,
    mean1 = 0.0,
    effect_size = 0.0,
    sd1 = 1,
    sd2 = 1,
    dist = "chi_square",
    par = NULL 
  )
  
  modifyList(defaults, list(...))
}

# get normality test object function

fn_to_get_norm_obj <- function(data) {
  return(data)
}

# downstream test 1(t-test) function

fn_for_ds_test_1 <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- t.test(x = x_data, y = y_data)
  return(list(p.value = test_result$p.value))
}

# downstream test 2(Mann-Whitney U Test) function
fn_for_ds_test_2 <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- wilcox.test(x_data, y_data)
  return(list(p.value = test_result$p.value))
}

