# -------Nonparametric: Mann-Whitney U Test-----
fn_for_ds_test_2 <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- wilcox.test(x_data, y_data)
  return(list(p.value = test_result$p.value))
}