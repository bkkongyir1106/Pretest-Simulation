#' Generate random samples from many named distributions
generate_samples <- function(n, dist) {
  # Normalize input and validate against the allowed set.
  # match.arg() will error with a helpful message if an unsupported name is given.
  dist <- match.arg(
    tolower(dist),
    c("normal", "normal_05", "normal_100","normal_50","normal_25","normal_15","normal_5","normal_30",
      "chi_square","chi_sq7","gamma","exponential", "t","t_5", "f", "uniform","uniform_100","uniform_10",
      "laplace","weibull","lognormal","contaminated","pareto","beta","cauchy","spike","extremeskew",
      "binomial","poisson","geometric","negative_binomial","inverse_gamma","inverse_gaussian","bimodal","triangular","logistic","gumbel","rayleigh")
  )
  
  # Each case returns a vector of length n.
  switch(dist,
         # --- Gaussian family (different locations/scales) ---
         normal              = rnorm(n, mean = 0,  sd = 1),
         normal_100          = rnorm(n, mean = 100, sd = 25),
         normal_50           = rnorm(n, mean = 50,  sd = 15),
         normal_25           = rnorm(n, mean = 25,  sd = 10),
         normal_15           = rnorm(n, mean = 15,  sd = 8),
         normal_5            = rnorm(n, mean = 5,   sd = 2),
         normal_30            = rnorm(n, mean = 30,   sd = 5),
         normal_05            = rnorm(n, mean = 10,   sd = 5),
         
         # --- Continuous, non-Gaussian ---
         chi_square          = rchisq(n, df = 3), 
         chi_sq7          = rchisq(n, df = 7),  
         gamma               = rgamma(n, shape = 2, rate = 1),        
         exponential         = rexp(n, rate = 1),                     
         
         # --- t distributions (heavy tails) ---
         t                   = rt(n, df = 3),
         t_5                 = rt(n, df = 5),
         
         # f distribution
         f = rf(n, df1 = 6, df2 = 15),
         
         # --- Uniforms over different ranges ---
         uniform             = runif(n, min = 0,  max = 1),
         uniform_100         = runif(n, min = 20, max = 100),
         uniform_10          = runif(n, min = 10, max = 50),
         
         # --- Symmetric, heavy-tailed / log-scale ---
         laplace             = LaplacesDemon::rlaplace(n, 0, 4),  
         weibull             = rweibull(n, shape = 2, scale = 1),
         lognormal           = rlnorm(n, meanlog = 0, sdlog = 1),
         
         # Mixture with inflated variance for 25% of draws (contamination)
         contaminated = {
           br <- rbinom(n, size = 1, prob = 0.75)       
           sd_br <- sqrt(1 + br * 24)                    
           rnorm(n, mean = 0, sd = sd_br)
         },
         
         pareto              = VGAM::rpareto(n = n, scale = 1, shape = 2), 
         beta                = rbeta(n, shape1 = 2, shape2 = 5),           
         cauchy              = rcauchy(n, location = 0, scale = 1),        
         
         # Point-mass “spike” plus normal slab
         spike = {
           spike_loc  <- 0
           spike_prop <- 0.1
           k <- floor(n * spike_prop)                    
           c(rep(spike_loc, k), rnorm(n - k))
         },
         
         # Extreme skew: cube of Exp() with random sign
         extremeskew = {
           y <- rexp(n)
           y^3 * sample(c(-1, 1), n, replace = TRUE)
         },
         
         # --- Discrete distributions ---
         binomial            = rbinom(n, size = 10, prob = 0.5),
         poisson             = rpois(n, lambda = 3),
         geometric           = rgeom(n, prob = 0.3),                   
         negative_binomial   = rnbinom(n, size = 5, prob = 0.5),
         
         # --- Inverse-family / special ---
         inverse_gamma       = LaplacesDemon::rinvgamma(n, shape = 3, scale = 2),
         inverse_gaussian    = statmod::rinvgauss(n, mean = 1, shape = 1),
         
         # --- Multimodal / triangular / logistic / extremes / Rayleigh ---
         bimodal             = c(rnorm(n %/% 2, mean = -2, sd = 1),
                                 rnorm(n - n %/% 2, mean =  2, sd = 1)),
         triangular          = extraDistr::rtriang(n, a = 0, b = 1, c = 0.5),
         logistic            = rlogis(n, location = 0, scale = 1),
         gumbel              = evd::rgumbel(n, loc = 0, scale = 1),
         rayleigh            = extraDistr::rrayleigh(n, sigma = 1)
  )
}
