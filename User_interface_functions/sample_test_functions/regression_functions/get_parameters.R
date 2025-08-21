get_parameters <- function(n, effect_size = NULL, ...) {
  defaults <- list(
    n = n,
    beta0 = 0,
    beta1 = 0,  # Default no effect
    x_dist = "Exponential",
    error_sd = 1,
    dist = "Normal"
  )
  
  args <- list(...)
  
  # If effect_size provided, use it for beta1
  if (!is.null(effect_size)) {
    args$beta1 <- effect_size
  }
  
  modifyList(defaults, args)
}