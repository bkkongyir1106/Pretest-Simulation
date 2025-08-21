#------- Kruskal-Wallis Test(Nonparameteric)----
fn_for_ds_test_2 <- function(data) {
  test_result <- kruskal.test(value ~ group, data = data)
  return(list(p.value = test_result$p.value))
}