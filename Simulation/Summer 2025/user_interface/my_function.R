# two-sample t-test
two_ind_t_test <- function(x, y = NULL, mu = 0, paired = FALSE) {
  t.test(x = x, y = y, mu = mu, paired = paired)
}

two_dep_t_test <- function(x, y = NULL, mu = 0, paired = TRUE) {
  t.test(x = x, y = y, mu = mu, paired = paired)
}

one_dep_t_test <- function(x, mu = 0) {
  t.test(x = x, mu = mu)
}
# linear regression
user_lm <- function(formula, data) {
  model <- lm(formula, data = data)
  return(summary(model))
}

# Mann Whithney U test
mw_u_test <- function(){
  wilcox.test(x, y)
}