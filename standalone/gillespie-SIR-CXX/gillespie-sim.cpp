/* ################################################################################
 *     _______ ____                _           _____ ____  _   __
 *    / ____(_) / /__  _________  (_)__       / ___// __ \/ | / /
 *   / / __/ / / / _ \/ ___/ __ \/ / _ \______\__ \/ /_/ /  |/ / 
 *  / /_/ / / / /  __(__  ) /_/ / /  __/_____/__/ / ____/ /|  /  
 *  \____/_/_/_/\___/____/ .___/_/\___/     /____/_/   /_/ |_/   
 *                      /_/      
 *                      
 *  Gillespie's direct method for a stochastic petri net (SPN) model
 *  Sean Wu (slwu89@berkeley.edu)
 *  May 2019
 * 
################################################################################ */


/* Rcpp header file */
#include <Rcpp.h>

/* algorithm for std::accumulate, iomanip for std::setw, string and sstream to write the message to break if h0 = 0 */
#include <algorithm>
#include <iomanip>
#include <string>
#include <sstream>

// [[Rcpp::export]]
Rcpp::NumericMatrix gillespie_CXX(const Rcpp::IntegerVector& M0,
                                 const double tmax,
                                 const Rcpp::List& pars,
                                 const Rcpp::Function& haz,
                                 const size_t info = 100,
                                 const size_t prealloc = 1e4
){
  
  /* set up local variables */
  double time = 0.;
  Rcpp::IntegerMatrix S = Rcpp::as<Rcpp::IntegerMatrix>(pars["S"]);
  Rcpp::IntegerVector M = M0;
  size_t u = S.nrow();
  size_t v = S.ncol();
  
  /* set up output */
  Rcpp::CharacterVector colnames = Rcpp::as<Rcpp::CharacterVector>(pars["P"]);
  colnames.push_front("time");
  Rcpp::NumericMatrix out(prealloc,u+1);
  std::fill(out.begin(),out.end(),R_NaN);
  
  /* record initial marking */
  size_t i = 0;
  out(i,0) = time;
  for(size_t j=0; j<M.size(); j++){
    out(i,j+1) = M[j];
  }
  
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
      std::stringstream msg;
      msg << " --- all events have probability 0, breaking from simulation at time: " << std::setw(4) << time << " --- \n";
      Rcpp::stop(msg.str());
    }
    
    /* Gillespie's direct method: P(what | when) * P(when) */
    time += R::rexp(1./h0);
    size_t j = Rcpp::as<size_t>(Rcpp::sample(v,1,false,h,false));
    
    /* update marking and output */
    M += S.column(j);
    i += 1;
    if(i >= out.nrow()){
      // make the new output
      Rcpp::NumericMatrix newout(prealloc*2,u+1);
      Rcpp::colnames(newout) = colnames;
      std::fill(newout.begin(),newout.end(),R_NaN);
      // copy the existing data into it
      for(size_t j=0; j<out.nrow(); j++){
        for(size_t k=0; j<out.ncol(); k++){
          newout(j,k) = out(j,k);
        }
      }
      out = newout;
    }
    
    out(i,0) = time;
    for(size_t j=0; j<M.size(); j++){
      out(i,j+1) = M[j];
    }
    
    if(i % info == 0){
      Rcpp::Rcout << " --- simulation at time: " << std::setw(4) << time << " --- \n";
    }
    
  }
  
  // only return the portion of the matrix with actual values
  Rcpp::NumericVector rowsums = Rcpp::rowSums(out);
  bool hasnan = Rcpp::any(Rcpp::is_nan(rowsums));
  if(hasnan){
    rowsums = Rcpp::na_omit(rowsums);
    out = out(Rcpp::Range(0,rowsums.size()-1), Rcpp::_ );
    Rcpp::colnames(out) = colnames;
    return out;
  } else {
    Rcpp::colnames(out) = colnames;
    return out;
  }
};