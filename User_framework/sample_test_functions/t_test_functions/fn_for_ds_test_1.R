# --------Parametric: Two-sample t-test --------
fn_for_ds_test_1 <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- t.test(x = x_data, y = y_data)
  return(list(p.value = test_result$p.value))
}
