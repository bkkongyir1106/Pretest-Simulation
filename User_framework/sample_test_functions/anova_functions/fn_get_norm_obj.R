#--------------- ANOVA residuals ----------------
fn_get_norm_obj <- function(data){
  return(residuals(aov(formula = value ~ group, data = data)))
}
