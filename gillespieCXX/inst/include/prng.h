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
 
// the question of why you would want to do this is nicely explained here:
// https://stackoverflow.com/questions/333889/why-have-header-files-and-cpp-files

// the reason you need header guards is explained surprisingly well at wikipedia: https://en.wikipedia.org/wiki/Include_guard
// and also here: https://stackoverflow.com/questions/2979384/purpose-of-header-guards
 
/* header guards */
#ifndef GILLESPIE_PRNG
#define GILLESPIE_PRNG

#include <Rcpp.h>
#include <random>
 
// a class that uses STL's Mersenne Twister (also what R uses, FWIW) to get "random" numbers ...
class prng {
public:
 
 /* constructor & destructor */
 prng(const uint_least32_t seed);
 ~prng(){};
 
 /* delete copy constructor/assignment operator, default move constructor/assignment operator */
 prng(const prng&) = delete;
 prng& operator=(const prng&) = delete;
 prng(prng&&) = default;
 prng& operator=(prng&&) = default;
 
 /* continuous random univariate sampling */
 double                                 get_runif();
 double                                 get_rexp(const double rate);
 
 /* discrete random multivariate sampling */
 size_t                                 get_rcategorical(const Rcpp::NumericVector& prob);
 
private:
 std::mt19937                            rng;
 std::uniform_real_distribution<double>  runif;
};

 
#endif