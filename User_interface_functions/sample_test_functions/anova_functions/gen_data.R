gen_data <- function(
    n = n,  
    effect_size = c(0, 0.2, 0.4, 0.5, 0.8),  
    sd = 1,            
    dist = "Normal"     
) {
  k <- length(effect_size)
  group_labels <- LETTERS[1:k]
  
  values <- unlist(lapply(effect_size, function(m) {
    m + sd * generate_data(n, dist)
  }))
  
  data.frame(
    group = rep(group_labels, each = n),
    value = values
  )
}