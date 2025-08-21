gen_data <- function(
    n1 = n,         
    n2 = n,         
    mean1 = 0.0,       
    effect_size = 0.0,     
    sd1 = 1,         
    sd2 = 1,         
    dist = "chi_square",
    par = NULL
) {
  
  group1 <- mean1 + sd1 * generate_data(n1, dist, par = par)
  group2 <- effect_size + sd2 * generate_data(n2, dist, par = par)
  
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
}

