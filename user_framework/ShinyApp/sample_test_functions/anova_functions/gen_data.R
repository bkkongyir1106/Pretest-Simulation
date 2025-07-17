gen_data <- function(
    n_per_group = 20,  
    means = c(0, 0.2, 0.4, 0.5, 0.8),  
    sd = 1,            
    dist = "LogNormal"     
) {
  k <- length(means)
  group_labels <- LETTERS[1:k]
  
  # Generate data for each group
  values <- unlist(lapply(means, function(m) {
    m + sd * generate_data(n_per_group, dist)
  }))
  
  return(data.frame(
    group = factor(rep(group_labels, each = n_per_group)),
    value = values
  ))
}
