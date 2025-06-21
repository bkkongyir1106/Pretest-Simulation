
# --------extract groups for normality test -------
fn_get_norm_obj <- function(data) {
  split(data$value, data$group)
}
