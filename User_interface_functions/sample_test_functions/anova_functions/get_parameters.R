get_parameters <- function(n, effect_size = NULL, ...) {
  defaults <- list(
    n = n,
    effect_size = c(0.0, 0.0, 0.0, 0.0, 0.0),
    sd = 1,
    dist = "Lognormal"
  )
  args <- list(...)
  if (!missing(effect_size)) {
    args$effect_size <- effect_size
  }
  modifyList(defaults, args) # Merge defaults with provided args
}