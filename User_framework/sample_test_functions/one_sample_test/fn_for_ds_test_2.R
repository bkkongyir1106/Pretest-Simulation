fn_for_ds_test_2 <- function(data) {
  test_result <- wilcox.test(data, mu = 0)
  return(list(p.value = test_result$p.value))
}