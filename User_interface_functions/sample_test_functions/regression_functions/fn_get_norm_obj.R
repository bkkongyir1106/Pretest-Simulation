fn_get_norm_obj <- function(data){
  model <- lm(y ~ x , data = data)
  return(residuals(model))
}

