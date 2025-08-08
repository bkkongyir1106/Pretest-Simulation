gen_data <- function(
    n1 = n,         
    n2 = n,         
    mean1 = 0.0,       
    means = 0.0,     
    sd1 = 1,         
    sd2 = 1,         
    dist = "LogNormal"  
) {
  
  group1 <- mean1 + sd1 * generate_data(n1, dist)
  group2 <- means + sd2 * generate_data(n2, dist)
  
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
}

