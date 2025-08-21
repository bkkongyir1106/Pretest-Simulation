gen_data <- function(
    n = 20,         
    effect_size = 0.5,       
    sd = 1,         
    dist = "Normal" ,
    par = NULL
) {
  
x <- effect_size + sd * generate_data(n, dist)
 
  return(x)
}

