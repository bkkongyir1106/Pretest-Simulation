  rm(list = ls())
  # setwd("/home/kongyir/spring2024/power")
  # source("/home/kongyir/spring2024/User_defined_functions.R")
  # source("/home/kongyir/spring2024/utility.R")
  setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research.05.10.2024/Power")
  source("~/Library/Mobile Documents/com~apple~CloudDocs/Research.05.10.2024/functions/User_defined_functions.R")
  source("~/Library/Mobile Documents/com~apple~CloudDocs/Research.05.10.2024/functions/utility.R")
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parallel process setup %%%%%%%%%%%%%%%%%%%----
  {
    par_set <- function(cores_reserve = 2) 
    {
      cores = parallel::detectCores()
      cores_use <- cores - cores_reserve
      if(Sys.info()["sysname"] == "Windows"){
        cl <- parallel::makeCluster(cores_use) # make a socket cluster
        doParallel::registerDoParallel(cl)     # for both Windows & Unix-like
        
      }else{
        
        cl <- snow::makeSOCKcluster(cores_use) # make a socket cluster
        doSNOW::registerDoSNOW(cl)           # for Unix or mac
      }
      foreach::getDoParWorkers()
      return(cl)
    }
    close_cluster <- function(cl) {
      parallel::stopCluster(cl) # close the cluster
    }
  }
  
  {## Set up
    N <- 1e3; threshold <- 0.05
    #dist_sum <- c("Standard Normal", "Uniform", "t" , "Contaminated", "Exponential", "Laplace", "Chi-Square", "Gamma", "Weibull", "LogNormal", "Pareto")
    dist_sum = "Pareto"
    nvec <- c(5, 10, 20, 25, 30) #, 40) 
    sig_level <- c( 0.05)
  }
  # Parallelized simulation setup
  {
    my_cl <- par_set(cores_reserve = 2)
    ntasks <- length(nvec) *  length(sig_level) * length(dist_sum) 
    pb <- txtProgressBar(max=ntasks, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  }
  # %%%%%%%%%%%%%%%%%%%%% SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%----
  system.time({
    sim_out <- foreach(n = nvec,
                       .packages = c("LaplacesDemon", "VGAM"),
                       .options.snow=opts) %:%
      foreach(dist = dist_sum) %:%
      foreach(alpha = sig_level) %dopar%
      {
        set.seed(1234)
        TotalSim.passed.SW.test = 0 ; TotalSim = 0 ; pval = c() 
        while (TotalSim.passed.SW.test < N) {
          x <- generate_data(n, dist) 
          y <- generate_data(n, dist) 
          TotalSim <- TotalSim + 1 
          if(shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha){
            TotalSim.passed.SW.test <- TotalSim.passed.SW.test + 1
            pval[TotalSim.passed.SW.test] <- t.test(x, y)$p.value
          }
        }
        error_t_test      <- mean(pval < alpha)
        prob.non.sig.SW.test <-   TotalSim.passed.SW.test/TotalSim
        Results <- list(
          error_t_test = error_t_test,
          prob.non.sig.SW.test = prob.non.sig.SW.test
        )
        return(Results)
      }
  })
  close_cluster(my_cl)        
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%% Output Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%-----
  powervec <- numeric(length(nvec) * length(dist_sum) * length(sig_level))
  TypeI.Error.T.Test <- prob.non.sig.SW.test <- array(powervec,dim = c(length(nvec), length(dist_sum), 
                            length(sig_level)), dimnames = list(nvec, dist_sum, sig_level))
  for (t in seq_along(nvec)) {
    for (j in seq_along(dist_sum)) {
      for (i in seq_along(sig_level)) {
        TypeI.Error.T.Test[t, j, i] <- (sim_out[[t]][[j]][[i]]$error_t_test)
        prob.non.sig.SW.test[t, j, i] <- (sim_out[[t]][[j]][[i]]$prob.non.sig.SW.test)
      }
    }
  }
  Inflation.TypeI.Eerror <- TypeI.Error.T.Test - threshold
  TypeI.Error.T.Test
  Inflation.TypeI.Eerror
  prob.non.sig.SW.test
  
  #save results
  #save.image(paste0("TwoSamplesExpectedInflation.TypeI.Errror20240608*",".RData"))

  save(nvec, TypeI.Error.T.Test, Inflation.TypeI.Eerror, prob.non.sig.SW.test, file = "Pareto_TwoSamplesExpectedInflation.TypeI.ErrorOR20240608.RData")
