/* ################################################################################
 *     _______ ____                _           _____ ____  _   __
 *    / ____(_) / /__  _________  (_)__       / ___// __ \/ | / /
 *   / / __/ / / / _ \/ ___/ __ \/ / _ \______\__ \/ /_/ /  |/ / 
 *  / /_/ / / / /  __(__  ) /_/ / /  __/_____/__/ / ____/ /|  /  
 *  \____/_/_/_/\___/____/ .___/_/\___/     /____/_/   /_/ |_/   
 *                      /_/      
 *                      
 *  Gillespie's direct method for a stochastic petri net (SPN) model
 *  a PRNG class
 *  Sean Wu (slwu89@berkeley.edu)
 *  May 2019
 * 
################################################################################ */
 
// include the header file so we can define how the class does what we have promised the compiler it can do!
#include "prng.h"
 

/* ################################################################################
* constructor
################################################################################ */

prng::prng(const uint_least32_t seed) : rng(seed){
  runif = std::uniform_real_distribution<double>(0,1);
};


/* ################################################################################
* continuous random univariate sampling
################################################################################ */

double prng::get_runif(){
  return runif(rng);
};

double prng::get_rexp(const double rate){
  std::exponential_distribution<double>rexp(rate);
  return rexp(rng);
};


/* ################################################################################
 * discrete random multivariate sampling
################################################################################ */
 
size_t prng::get_rcategorical(const Rcpp::NumericVector& prob){
 std::discrete_distribution<int>rcategorical(prob.begin(),prob.end());
 return rcategorical(rng);
};
