gen_data <- function(
    n = 20,         
    mean = 0.5,       
    sd = 1,         
    dist = "Exponential" ,
    par = NULL
) {
  
x <- mean + sd * generate_data(n, dist)
 
  return(x)
}

