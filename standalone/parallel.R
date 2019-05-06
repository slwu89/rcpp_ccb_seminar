################################################################################
#      _____ ________
#     / ___//  _/ __ \
#     \__ \ / // /_/ /
#    ___/ // // _, _/
#   /____/___/_/ |_|
# 
#   How to do this in parallel with foreach
#   Sean Wu (slwu89@berkeley.edu)
#   May 2019
# 
################################################################################


################################################################################
#   load packages and source files
################################################################################

rm(list=ls());gc()

library(here)
library(Rcpp)

library(parallel)
library(foreach)
library(doSNOW)
library(abind)

library(ggplot2)
library(reshape2)
library(matrixStats)


################################################################################
#   summarize a lot of Monte Carlo runs
################################################################################

# discretise continuous time output of the CTMC
# adapted from the "smfsb" package
discretise <- function(out, start =0, end = NaN, dt=1){
  
  events <- length(out[,"time"])
  
  if(is.nan(end)){
    end <- unname(out[events,"time"])
  }
  if(!is.nan(end) & end < out[events,"time"]){
    out <- out[out[,"time"] <= end,]
    events <- length(out[,"time"])
  }
  
  len <- (end-start)%/%dt+1
  x <- matrix(0,nrow=len,ncol=ncol(out),dimnames = dimnames(out))
  x[,"time"] <- (0:(len-1))*dt
  
  # for trajectories that die off immediately
  if(nrow(out) < 2){
    x[i,-1] <- out[1,-1]
    return(x)
  }
  
  target <- 0
  j <- 1
  
  for (i in 1:events) {
    while (out[i,"time"] >= target) {
      x[j,-1] <- out[i,-1]
      j <- j+1
      target <- target+dt
    }
  }
  
  return(x)
}


# make the SPN
SPN_T <- c("inf","rec")
SPN_P <- c("S","I","R")

Pre <- matrix(data = c(
  1,1,0,
  0,1,0
),nrow = length(SPN_T),ncol = length(SPN_P),byrow = TRUE,
dimnames = list(SPN_T,SPN_P))

Post <- matrix(data = c(
  0,2,0,
  0,0,1
),nrow = length(SPN_T),ncol = length(SPN_P),byrow = TRUE,
dimnames = list(SPN_T,SPN_P))

A <- Post - Pre
S <- t(A)

# function to return vector of hazards
hazard <- function(t,state,pars,...){
  
  with(pars,{
    # browser()
    h <- rep(0,length(T))
    h <- setNames(h,T)
    
    enabled <- Pre %*% state
    
    # infection dynamics
    if(enabled["inf",] > 0){
      h["inf"] <- state["S"]*state["I"]*beta
    }
    
    if(enabled["rec",] > 0){
      h["rec"] <- state["I"]*gamma
    }
    
    return(h)  
  })
  
}

M0 <- setNames(c(1e3,1,0),SPN_P)

theta <- list("beta"=.0001,"gamma"=1/50)
pars <- list(Pre=Pre,Post=Post,S=S,T=SPN_T,P=SPN_P)
pars <- c(pars,theta)

nrep <- 1e3
tmax <- 250

# set up cluster and source the file on each core
# if you don't do this you will likely get a strange error about a NULL value passed as symbol address, because reasons
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)
clusterEvalQ(cl,{
  Rcpp::sourceCpp(here::here("gillespie-SIR-CXX/gillespie-sim.cpp"))
})

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

# combine each matrix into a slice of a 3d array
acomb <- function(...) abind(..., along=3)

SIR_ensemble <- foreach(i = 1:nrep,.combine = acomb,.options.snow=opts) %dopar% {
  
  SIR <- gillespie_CXX(M0 = M0,tmax = tmax,pars = pars,haz = hazard,info = 1e6)
  discretise(out = SIR,start = 0,end = tmax-1,dt = 1)
  
}

close(pb)
stopCluster(cl);rm(cl);gc()

SIR_means <- data.frame(
  time = 0:(tmax-1),
  S = rowMeans(SIR_ensemble[,"S",]),
  I = rowMeans(SIR_ensemble[,"I",]),
  R = rowMeans(SIR_ensemble[,"R",])
)

SIR_means <- melt(SIR_means,id.vars = "time")

SIR_sds <- data.frame(
  time = 0:(tmax-1),
  S = rowSds(SIR_ensemble[,"S",]),
  I = rowSds(SIR_ensemble[,"I",]),
  R = rowSds(SIR_ensemble[,"R",])
)

SIR_sds <- melt(SIR_sds,id.vars = "time")

SIR_out <- merge(x = SIR_means,y = SIR_sds,by = c("time","variable"),suffixes = c("mean","sd"),sort = F)

ggplot(data = SIR_out) +
  geom_line(aes(x=time,y=valuemean,color=variable)) +
  geom_ribbon(aes(x=time,ymin=valuemean-valuesd,ymax=valuemean+valuesd,fill=variable),alpha=0.25) +
  theme_bw()

# low <- 0.2
# high <- 0.8
# 
# SIR_quant_low <- data.frame(
#   time = 0:(tmax-1),
#   S = apply(SIR_ensemble[,"S",],MARGIN = 1,function(x){quantile(x,probs=(low))} ),
#   I = apply(SIR_ensemble[,"I",],MARGIN = 1,function(x){quantile(x,probs=(low))} ),
#   R = apply(SIR_ensemble[,"R",],MARGIN = 1,function(x){quantile(x,probs=(low))} )
# )
# SIR_quant_low <- melt(SIR_quant_low,id.vars = "time")
# 
# SIR_quant_high <- data.frame(
#   time = 0:(tmax-1),
#   S = apply(SIR_ensemble[,"S",],MARGIN = 1,function(x){quantile(x,probs=(high))} ),
#   I = apply(SIR_ensemble[,"I",],MARGIN = 1,function(x){quantile(x,probs=(high))} ),
#   R = apply(SIR_ensemble[,"R",],MARGIN = 1,function(x){quantile(x,probs=(high))} )
# )
# SIR_quant_high <- melt(SIR_quant_high,id.vars = "time")
# 
# SIR_quant <- merge(SIR_quant_low,SIR_quant_high,by = "time",suffixes = c("low","high"))
# 
# ggplot() +
#   geom_line(aes(x=time,y=value,color=variable), data =SIR_means) +
#   geom_ribbon(aes(x=time,ymin=valuelow,ymax=valuehigh,fill=variablelow,group=variablelow),alpha=0.2,data = SIR_quant) +
#   theme_bw()
