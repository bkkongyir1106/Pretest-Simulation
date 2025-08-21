get_parameters <- function(n, ...) {
  defaults <- list(
    n1 = n,
    n2 = n,
    mean1 = 0.0,
    effect_size = 0.0,
    sd1 = 1,
    sd2 = 1,
    dist = "chi_square",
    par = NULL
  )
  
  modifyList(defaults, list(...))
}
