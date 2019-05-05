################################################################################
#      _____ ________        ____ 
#     / ___//  _/ __ \      / __ \
#     \__ \ / // /_/ /_____/ /_/ /
#    ___/ // // _, _/_____/ _, _/ 
#   /____/___/_/ |_|     /_/ |_|  
# 
#   The stochastic SIR model example
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

source(here::here("gillespie-SIR-R/gillespie-sim.R"))
sourceCpp(here::here("gillespie-SIR-CXX/gillespie-sim.cpp"))


################################################################################
#   create the SPN model
################################################################################

# initial marking of net
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
    
    # if(enabled["wane",] > 0){
    #   h["wane"] <- state["R"]*delta
    # }
    
    return(h)  
  })
  
}


################################################################################
#   run the stochastic SIR model
################################################################################

# initial marking of SPN
M0 <- setNames(c(1e3,1,0),SPN_P)

theta <- list("beta"=.0001,"gamma"=1/50,"delta"=1/365)
pars <- list(Pre=Pre,Post=Post,S=S,T=SPN_T,P=SPN_P)
pars <- c(pars,theta)

set.seed(42)
SIR <- gillespie_CXX(M0 = M0,tmax = 250,pars = pars,haz = hazard,info = 100)
SIR <- gillespie_R(M0 = M0,tmax = 250,pars = pars,haz = hazard,info = 100)

# plot trajectory
ymax <- max(SIR[,2:ncol(SIR)]) + 10
SIRcolor <-   c(S="steelblue",I="firebrick3",R="darkorchid3")
plot(x = SIR[,"time"],y = SIR[,"S"],type = "l",col = SIRcolor["S"],lwd = 1.85,
     ylim = c(0,ymax),xlab = "Time (days)",ylab = "Count",main = "Stochastic SIR")
lines(x = SIR[,"time"],y = SIR[,"I"],col = SIRcolor["I"],lwd = 1.85)
lines(x = SIR[,"time"],y = SIR[,"R"],col = SIRcolor["R"],lwd = 1.85)

# just plot the infecteds
# plot(x = SIR[,"time"],y = SIR[,"I"],type = "l",col = SIRcolor["I"],lwd = 1.5,
#      xlab = "Time (days)",ylab = "Infecteds")


################################################################################
#   run the version using header files
################################################################################

sourceCpp(here::here("gillespie-SIR-CXX-headers/gillespie-sim-headers.cpp"))

# please don't use this to get seeds for real applications
seed <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)

SIR <- gillespie_CXX_h(M0 = M0,tmax = 250,pars = pars,haz = hazard,seed = seed,info = 100)

ymax <- max(SIR[,2:ncol(SIR)]) + 10
SIRcolor <-   c(S="steelblue",I="firebrick3",R="darkorchid3")
plot(x = SIR[,"time"],y = SIR[,"S"],type = "l",col = SIRcolor["S"],lwd = 1.85,
     ylim = c(0,ymax),xlab = "Time (days)",ylab = "Count",main = "Stochastic SIR")
lines(x = SIR[,"time"],y = SIR[,"I"],col = SIRcolor["I"],lwd = 1.85)
lines(x = SIR[,"time"],y = SIR[,"R"],col = SIRcolor["R"],lwd = 1.85)