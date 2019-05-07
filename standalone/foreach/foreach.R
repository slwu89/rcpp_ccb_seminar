Rcpp::sourceCpp(here::here('foreach/timestwo.cpp'))

timesTwo(5)

library(parallel)
library(foreach)
library(doSNOW)

cl <- makeSOCKcluster(2)
registerDoSNOW(cl)

# rm(timesTwo)
# 
# clusterEvalQ(cl,{
#   Rcpp::sourceCpp(here::here('foreach/timestwo.cpp'))
# })


foreach(i = 1:10,.combine = "c") %dopar% {
  timesTwo(i)
}

stopCluster(cl)
rm(cl);gc()
