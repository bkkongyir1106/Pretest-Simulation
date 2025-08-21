get_parameters <- function(n, ...) {
  defaults <- list(
    n_per_group = n,
    means = c(0.0, 0.0, 0.0, 0.0, 0.0),
    sd = 1,
    dist = "Normal"
  )
  modifyList(defaults, list(...))
}