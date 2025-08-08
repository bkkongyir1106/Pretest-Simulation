gen_data <- function(
    n_per_group = n,  
    means = c(0, 0.2, 0.4, 0.5, 0.8),  
    sd = 1,            
    dist = "Normal"     
) {
  k <- length(means)
  group_labels <- LETTERS[1:k]
  
  values <- unlist(lapply(means, function(m) {
    m + sd * generate_data(n_per_group, dist)
  }))
  
  data.frame(
    group = rep(group_labels, each = n_per_group),
    value = values
  )
}