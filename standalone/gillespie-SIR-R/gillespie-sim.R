################################################################################
#      _______ ____                _           _____ ____  _   __
#     / ____(_) / /__  _________  (_)__       / ___// __ \/ | / /
#    / / __/ / / / _ \/ ___/ __ \/ / _ \______\__ \/ /_/ /  |/ / 
#   / /_/ / / / /  __(__  ) /_/ / /  __/_____/__/ / ____/ /|  /  
#   \____/_/_/_/\___/____/ .___/_/\___/     /____/_/   /_/ |_/   
#                       /_/                                      
# 
#   Gillespie's direct method for a stochastic petri net (SPN) model
#   Sean Wu (slwu89@berkeley.edu)
#   May 2019
# 
################################################################################

# Gillespie simulator for a SPN
gillespie_R <- function(M0,tmax,pars,haz,info=100,prealloc=1e4){
  
  P <- pars$P
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