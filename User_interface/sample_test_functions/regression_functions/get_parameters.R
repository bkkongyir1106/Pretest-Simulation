get_parameters <- function(
    n        = n,
    beta0     = 0,
    means     = 0.0,           
    x_dist    = "Exponential",
    error_sd  = 1,
    dist= "Normal",
    ...
) {
  defaults <- list(
    n         = n,
    beta0     = beta0,
    means     = means,
    x_dist    = x_dist,
    error_sd  = error_sd,
    dist      = dist
  )
  modifyList(defaults, list(...))
}


