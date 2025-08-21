gen_data <- function(
    n,           
    beta0 = 0,        
    beta1 = effect_size, 
    x_dist = "Exponential",
    error_sd = 1,     
    dist = "Normal"  
) {
  # predictor
  x <- generate_data(n, dist = x_dist)
  # error
  error <- error_sd * generate_data(n, dist)
  # independent variable
  y <- beta0 + beta1 * x + error
  return(data.frame(x = x, y = y))
}

