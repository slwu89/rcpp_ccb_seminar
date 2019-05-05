#include <Rcpp.h>
#include <algorithm>
#include <iomanip>

// [[Rcpp::export]]
Rcpp::NumericMatrix gillespieCXX(const Rcpp::IntegerVector& M0,
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
  Rcpp::colnames(out) = colnames;
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
      Rcpp::Rcout << " --- all events have probability 0, breaking from simulation at time: " << std::setw(4) << time << " --- \n";
      Rcpp::stop("\n");
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
  
  // Rcpp::IntegerVector rowsums = Rcpp::rowSums(out);
  // Rcpp::LogicalVector row2r = !Rcpp::is_nan(rowsums);
  return out;
};