# -------- Define ANOVA parameter generator -----
get_parameters <- function(n) {
  list(
    n_per_group = n,
    means = c(0.0, 0.5, 0.8),
    sd = 1,
    dist = "Normal"
  )
}
