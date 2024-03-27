#%%%%%%%%%%%%%%%%%%%%% (1) UTILS %%%%%%%%%%%%%%%%%%%%%%%----
if(!require("pacman")) install.packages("pacman")
pacman::p_load(MASS, doMC, microbenchmark, psych, kernlab, tidyverse, haven, foreach, doParallel, doRNG,
               Matrix, latex2exp, doSNOW, nortest, dgof, moments, stats, tseries, LaplacesDemon)

## import utility files
path = "~/Documents/Testing Files/ExpectileGEE"                  # set path
file.sources = list.files(path,                                  # locate all .R files
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)
