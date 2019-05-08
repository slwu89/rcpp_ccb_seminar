// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gillespie_dm
Rcpp::DataFrame gillespie_dm(const Rcpp::IntegerVector& M0, const double tmax, const Rcpp::List& pars, const Rcpp::Function& haz, const uint_least32_t seed, const size_t info, const size_t prealloc);
RcppExport SEXP _gillespieCXX_gillespie_dm(SEXP M0SEXP, SEXP tmaxSEXP, SEXP parsSEXP, SEXP hazSEXP, SEXP seedSEXP, SEXP infoSEXP, SEXP preallocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type M0(M0SEXP);
    Rcpp::traits::input_parameter< const double >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type haz(hazSEXP);
    Rcpp::traits::input_parameter< const uint_least32_t >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const size_t >::type info(infoSEXP);
    Rcpp::traits::input_parameter< const size_t >::type prealloc(preallocSEXP);
    rcpp_result_gen = Rcpp::wrap(gillespie_dm(M0, tmax, pars, haz, seed, info, prealloc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gillespieCXX_gillespie_dm", (DL_FUNC) &_gillespieCXX_gillespie_dm, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_gillespieCXX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}