fn_for_ds_test_1<- function(data) {
  model <- lm(y ~ x, data = data)
  tidy_model <- broom::tidy(model)
  p_value <- tidy_model$p.value[tidy_model$term == "x"]
  return(list(p.value = p_value))
}
