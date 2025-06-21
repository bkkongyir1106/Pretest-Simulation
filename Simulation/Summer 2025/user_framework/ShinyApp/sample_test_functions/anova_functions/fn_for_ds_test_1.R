#------------- ANOVA Test ---------------------
fn_for_ds_test_1 <- function(data) {
  aov_model <- aov(value ~ group, data = data)
  p_value <- summary(aov_model)[[1]]$"Pr(>F)"[1]
  return(list(p.value = p_value))
}

