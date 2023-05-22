// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mfvb
Rcpp::List mfvb(arma::vec response, arma::mat fir_l, arma::mat sec_l, double a_tau, double b_tau, double a_pi, double b_pi, int max_iter, double max_tol);
RcppExport SEXP _GraphR_mfvb(SEXP responseSEXP, SEXP fir_lSEXP, SEXP sec_lSEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP a_piSEXP, SEXP b_piSEXP, SEXP max_iterSEXP, SEXP max_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type response(responseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fir_l(fir_lSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sec_l(sec_lSEXP);
    Rcpp::traits::input_parameter< double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< double >::type b_tau(b_tauSEXP);
    Rcpp::traits::input_parameter< double >::type a_pi(a_piSEXP);
    Rcpp::traits::input_parameter< double >::type b_pi(b_piSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type max_tol(max_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(mfvb(response, fir_l, sec_l, a_tau, b_tau, a_pi, b_pi, max_iter, max_tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GraphR_mfvb", (DL_FUNC) &_GraphR_mfvb, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_GraphR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
