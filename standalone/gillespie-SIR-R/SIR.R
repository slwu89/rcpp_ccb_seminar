################################################################################
#      _____ ________        ____ 
#     / ___//  _/ __ \      / __ \
#     \__ \ / // /_/ /_____/ /_/ /
#    ___/ // // _, _/_____/ _, _/ 
#   /____/___/_/ |_|     /_/ |_|  
# 
#   The stochastic SIR model example in R
#   Sean Wu (slwu89@berkeley.edu)
#   May 2019
# 
################################################################################

rm(list=ls());gc()

# initial marking of net
T <- c("inf","rec")
P <- c("S","I","R")

Pre <- matrix(data = c(
  1,1,0,
  0,1,0
),nrow = length(T),ncol = length(P),byrow = TRUE,
dimnames = list(T,P))

Post <- matrix(data = c(
  0,2,0,
  0,0,1
),nrow = length(T),ncol = length(P),byrow = TRUE,
dimnames = list(T,P))

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
      beta <- R0*gamma
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

# Gillespie simulator for a SPN
gillespie <- function(M0,tmax,pars,haz,info=100,prealloc=1e4){
  
  S <- pars$S
  time <- 0
  M <- M0
  u <- nrow(S)
  v <- ncol(S)
  
  out <- matrix(data = NaN,nrow = prealloc,ncol = u+1,dimnames = list(NULL,c("time",P)))
  i <- 1
  out[i,] <- c(time,M)
  
  while(time < tmax){
    
    # evaluate hazard function
    h <- haz(time,M,pars)
    
    # if all events have probability 0, break the loop
    if(sum(h) <= .Machine$double.eps){
      cat(" --- all events have probability 0, breaking from simulation at time: ",round(time,3)," --- \n")
      break
    }
    
    # Gillespie's direct method: P(what | when) * P(when)
    time <- time + rexp(n = 1,rate = sum(h))
    j <- sample.int(n = v,size = 1,prob = h)
    
    # update marking and output
    M <- M + S[,j]
    i <- i + 1
    if(i > nrow(out)){
      out <- rbind(out,matrix(data = NaN,nrow = prealloc,ncol = u+1,dimnames = list(NULL,c("time",P))))
    }
    out[i,] <- c(time,M)
    
    if(i%%info==0){
      cat(" --- simulation at time: ",round(time,3)," --- \n")
    }
  }
  
  # dont return rows that are just for memory pre-allocation
  return(
    out[!is.nan(rowSums(out)),]
  )
}

# run the stochastic epidemic model
M0 <- setNames(c(1e3,1,0),P)

theta <- list("R0"=0.005,"gamma"=1/50,"delta"=1/365)
pars <- list(Pre=Pre,Post=Post,S=S,T=T,P=P)
pars <- c(pars,theta)

SIR <- gillespieCXX(M0 = M0,tmax = 250,pars = pars,haz = hazard,info = 100)
SIR <- SIR[!is.na(rowSums(SIR)),]

SIR <- gillespie(M0 = M0,tmax = 250,pars = pars,haz = hazard,info = 100)

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