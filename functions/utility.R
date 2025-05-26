#%%%%%%%%%%%%%%%%%%%%% (1) UTILS %%%%%%%%%%%%%%%%%%%%%%%----
if(!require("pacman")) install.packages("pacman")
pacman::p_load(MASS, doMC, microbenchmark, psych, kernlab, tidyverse,
               haven, foreach, doParallel, doRNG, cramer,Matrix, latex2exp, 
               doSNOW, nortest, dgof, goftest, moments, stats, tseries, LaplacesDemon,
               VGAM, e1071, openxlsx, BayesFactor, boot, rstan, reshape2, ggplot2)


