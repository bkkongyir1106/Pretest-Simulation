gen_data <- function(distn, sample_size, ...) {
  return(do.call(distn, list(sample_size,...)))
}

not_normality_test <- function(sample_data, shapiro.alpha) {
  return(shapiro.test(sample_data)$p.value <= shapiro.alpha)
}

sim <- function(iter_passed_size, distn, mu, shapiro.alpha, t.alpha, sample_size, ...){
  count = 0  #No of times data is generated irrespective of passing normality
  reject_h0 = 0 # No fo times we reject the null hypothesis
  total_s.test_passed = iter_passed_size
  while(iter_passed_size > 0) {
    x <- gen_data(distn, sample_size, ...)
    count = count + 1
    if(not_normality_test(x,shapiro.alpha)) {
      if(wilcox.test(x, mu=mu)$p.value < t.alpha) {
        reject_h0 = reject_h0 + 1
      }
      iter_passed_size = iter_passed_size - 1
    }
  }
  print(c("No of Data Generated: " , count))
  print(c("No of times Normality passed: ", total_s.test_passed))
  print(c("T-Test Type 1 error rate(alpha: ", reject_h0 / total_s.test_passed))
}
sim(10000,"rexp",1,0.05, 0.05,10)
sim(10000,"runif",1/2,0.05, 0.05,10)
sim(10000,"rnorm",0,0.05, 0.05,10)
