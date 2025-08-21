fn_for_ds_test_1 <- function(data) {
  test_result <- t.test(data, mu = 0)
  return(list(p.value = test_result$p.value))
}
