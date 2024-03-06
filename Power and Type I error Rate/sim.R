
num.sim <- 1e3
n <- 10
df <- 5
deltavec <- c(0, 1 ,2, 3)
tstat <- tstaty <- matrix(NA, nrow = num.sim, ncol = length(deltavec))

set.seed(123)

for(i in 1:num.sim){
  for (k in 1:length(deltavec)){
    x <- (rchisq(n, df) - df)/sqrt(2*df) + deltavec[k]
    y <- rnorm(n) + deltavec[k]
    tstat[i, k] <- sqrt(n)*mean(x)/sd(x)
    tstaty[i, k] <- sqrt(n)*mean(y)/sd(y)
  }
}

plot(density(tstat[, 1]), typ = "n", xlim = c(-5, 15))
legend_text <- c()
for (k in 1:length(deltavec)){
  lines(density(tstaty[, k]), col = "red", lty = k)
  lines(density(tstat[, k]), col = "blue", lty = k)
  legend_text <- c(legend_text, col = k) 
}
legend('topright', legend = c("chi-squared", "normal"), col=c("red", "blue"))


