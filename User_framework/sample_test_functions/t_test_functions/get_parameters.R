get_parameters <- function(n,...) {
  defaults <- list(
    n1 = n,
    n2 = n,
    mean1 = 0.0,
    means = 0.0,
    sd1 = 1,
    sd2 = 1,
    dist = "Uniform"
  )
  
  modifyList(defaults, list(...))
}
