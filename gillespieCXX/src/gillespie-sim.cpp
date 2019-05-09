/* ################################################################################
 *     _______ ____                _           _____ ____  _   __
 *    / ____(_) / /__  _________  (_)__       / ___// __ \/ | / /
 *   / / __/ / / / _ \/ ___/ __ \/ / _ \______\__ \/ /_/ /  |/ / 
 *  / /_/ / / / /  __(__  ) /_/ / /  __/_____/__/ / ____/ /|  /  
 *  \____/_/_/_/\___/____/ .___/_/\___/     /____/_/   /_/ |_/   
 *                      /_/      
 *                      
 *  Gillespie's direct method for a stochastic petri net (SPN) model
 *  Using STL's random number generator in a custom class
 *  Sean Wu (slwu89@berkeley.edu)
 *  May 2019
 * 
################################################################################ */
 
 
/* Rcpp header file */
#include <Rcpp.h>

/* algorithm for std::accumulate, iomanip for std::setw, string and sstream to write the message to break if h0 = 0 */
#include <algorithm>
#include <iomanip>
#include <sstream>

/* for output */
#include <vector> 

/* for unique_ptr */
#include <memory>
 
/* the class for our random numbers */
#include "prng.h"

using prng_ptr = std::unique_ptr<prng>;
 
//' Gillespie's Direct Method in C++
//' 
//' Simulate a trajectory from a stochastic petri net model.
//' 
//' @param M0 the initial marking of the net
//' @param tmax the maximum simulation time
//' @param pars a named list of parameters to pass to \code{haz}
//' @param haz an R function that takes 3 arguments (time, state, pars) and returns a vector of hazards for each event
//' @param seed a positive integer seed for the random number generator (\url{http://www.cplusplus.com/reference/random/})
//' @param info print information on simulation every \code{info} events
//' @param prealloc how much memory to preallocate for output
//' 
//' @examples
//' \dontrun{
//' library(gillespieCXX)
//' library(ggplot2)
//' library(reshape2)
//' data(SIR)
//' out <- gillespie_dm(M0 = SIR$M0,tmax = 250,pars = SIR$pars,haz = SIR$hazard,seed = 12214L,info = 100)
//' out <- as.data.frame(out)
//' ggplot(data = melt(out,id.vars="time")) +
//'   geom_line(aes(x=time,y=value,color=variable)) +
//'   theme_bw()
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame gillespie_dm(const Rcpp::IntegerVector& M0,
                                 const double tmax,
                                 const Rcpp::List& pars,
                                 const Rcpp::Function& haz,
                                 const uint_least32_t seed,
                                 const size_t info = 100,
                                 const size_t prealloc = 1e4
){
 
 /* initialize prng */
 prng_ptr PRNG = std::make_unique<prng>(seed);
 
 /* set up local variables */
 double time = 0.;
 Rcpp::IntegerMatrix S = Rcpp::as<Rcpp::IntegerMatrix>(pars["S"]);
 Rcpp::IntegerVector M = M0;
 size_t u = S.nrow();
 
 /* set up output */
 std::vector<std::vector<double> > vecout(u+1);
 for(size_t i=0; i<vecout.size(); i++){
   vecout.at(i).reserve(prealloc);
 }
 
 /* record initial marking */
 vecout.at(0).emplace_back(time);
 for(size_t j=0; j<u; j++){
   vecout.at(j+1).emplace_back((double)M[j]);
 }
 
 /* used to schedule how often to check for aborts */
 size_t i = 0;
 
 /* main simulation loop */
 while(time < tmax){
   
   /* need to check if a user has pressed abort from R's REPL */
   if(i % 100 == 0){
     Rcpp::checkUserInterrupt();
   }
   
   /* evaluate hazard function */
   Rcpp::NumericVector h = haz(time,M,pars);
   
   /* if all events have probability 0, break the loop */
   double h0 = std::accumulate(h.begin(),h.end(),0.);
   if(h0 <= 2.22e-16){
     Rcpp::Rcout << " --- all events have probability 0, breaking from simulation at time: " << std::setw(4) << time << " --- \n";
     break;
   }
   
   /* Gillespie's direct method: P(what | when) * P(when) */
   time += PRNG->get_rexp(h0);
   size_t j = PRNG->get_rcategorical(h);
   
   /* update marking and output */
   M += S.column(j);
   i += 1;
   
   vecout.at(0).emplace_back(time);
   for(size_t j=0; j<u; j++){
     vecout.at(j+1).emplace_back(M[j]);
   }
   
   if(i % info == 0){
     Rcpp::Rcout << " --- simulation at time: " << std::setw(4) << time << " --- \n";
   }
   
 }
 
 /* return output */
 Rcpp::CharacterVector colnames = Rcpp::as<Rcpp::CharacterVector>(pars["P"]);
 colnames.push_front("time");
 
 Rcpp::DataFrame output;
 for(size_t j=0; j<vecout.size(); j++){
   output.push_back(Rcpp::wrap(vecout.at(j)));
 }
 output.attr("names") = colnames;
 return output;
};