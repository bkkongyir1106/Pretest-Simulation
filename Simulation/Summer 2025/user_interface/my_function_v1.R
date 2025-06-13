calculate_residuals <- function(model){
  return(resid(model))
}

calculate_samples <- function(data) {
  if (is.data.frame(data) || is.matrix(data)) {
    return(as.list(as.data.frame(data)))
  } else if (is.list(data)) {
    return(data)
  } else {
    return(list(data))  
  }
}

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
Mann_Whitney_U_test <- function(x, y = NULL, mu = 0){
  return(wilcox.test(x, y, mu = mu))
}

# ANOVA
anova_test <- function(data){
  return(residuals(aov(formula = value ~ group, data = all_data())))
}


# example data generation functions

data_gen <- function(){
  data <- data.frame(
             x = rnorm(10, mean = 0, sd = 1),
             y = rexp(10, rate = 1),
             z = rchisq(10, df = 3))
  return(data)
}

n <- 50
all_data <- function(){
  data <- list(a = generate_data(n = n, dist = "Gamma"),
       b = generate_data(n = n, dist = "Uniform"),
       c = generate_data(n = n, dist = "Normal"),
       d = generate_data(n = n, dist = "Normal"),
       x = generate_data(n = n, dist = "Normal"),
       y = generate_data(n = n, dist = "Exponential"),
       z = generate_data(n = n, dist = "LogNormal"))
  
  # Combine into a data frame for anova
  group <- factor(rep(c("a","b" ,"x", "y", "z"), each = n))
  value <- c(data$a, data$b, data$x, data$y, data$z)
   return(data.frame(group, value))
}



