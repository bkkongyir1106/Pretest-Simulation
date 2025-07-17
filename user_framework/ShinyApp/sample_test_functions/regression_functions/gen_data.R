gen_data <- function(
    n = 30,           
    beta0 = 0,        
    beta1 = 0.5, 
    x_dist = "Exponential",
    error_sd = 1,     
    error_dist = "Normal"  
) {
  # predictor
  x <- generate_data(n, dist = x_dist)
  # error
  error <- error_sd * generate_data(n, error_dist)
  # independent variable
  y <- beta0 + beta1 * x + error
  return(data.frame(x = x, y = y))
}

