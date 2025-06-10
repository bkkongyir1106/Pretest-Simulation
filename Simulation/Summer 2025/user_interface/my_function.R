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