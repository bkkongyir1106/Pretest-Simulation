gen_data <- function(
    n = n,           
    beta0 = 0,        
    means = 0.0, 
    x_dist = "Exponential",
    error_sd = 1,     
    dist = "Normal"  
) {
  # predictor
  x <- generate_data(n, dist = x_dist)
  # error
  error <- error_sd * generate_data(n, dist)
  # independent variable
  #beta1 <- means
  y <- beta0 + means * x + error
  return(data.frame(x = x, y = y))
}

