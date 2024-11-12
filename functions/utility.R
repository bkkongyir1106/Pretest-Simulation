#%%%%%%%%%%%%%%%%%%%%% (1) UTILS %%%%%%%%%%%%%%%%%%%%%%%----
if(!require("pacman")) install.packages("pacman")
pacman::p_load(MASS, doMC, microbenchmark, psych, kernlab, tidyverse,
               haven, foreach, doParallel, doRNG, cramer,Matrix, latex2exp, 
               doSNOW, nortest, dgof, goftest, moments, stats, tseries, LaplacesDemon,
               VGAM, e1071, openxlsx, BayesFactor)

## import utility files
path = "~/Documents/Testing Files/ExpectileGEE"                  # set path
file.sources = list.files(path,                                  # locate all .R files
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)


