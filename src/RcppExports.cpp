// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hungarian_cc
arma::umat hungarian_cc(arma::mat cost);
RcppExport SEXP _COMIX_hungarian_cc(SEXP costSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cost(costSEXP);
    rcpp_result_gen = Rcpp::wrap(hungarian_cc(cost));
    return rcpp_result_gen;
END_RCPP
}
// calib
Rcpp::List calib(arma::mat Y, arma::vec C, arma::mat Z, NumericVector mu_input, IntegerVector mu_dim, NumericVector mu0_input, IntegerVector mu0_dim, int ref);
RcppExport SEXP _COMIX_calib(SEXP YSEXP, SEXP CSEXP, SEXP ZSEXP, SEXP mu_inputSEXP, SEXP mu_dimSEXP, SEXP mu0_inputSEXP, SEXP mu0_dimSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_input(mu_inputSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mu_dim(mu_dimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0_input(mu0_inputSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mu0_dim(mu0_dimSEXP);
    Rcpp::traits::input_parameter< int >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(calib(Y, C, Z, mu_input, mu_dim, mu0_input, mu0_dim, ref));
    return rcpp_result_gen;
END_RCPP
}
// calibNoDist
Rcpp::List calibNoDist(arma::mat Y, arma::vec C, arma::mat Z, NumericVector mu_input, IntegerVector mu_dim, NumericVector mu0_input, IntegerVector mu0_dim, int ref);
RcppExport SEXP _COMIX_calibNoDist(SEXP YSEXP, SEXP CSEXP, SEXP ZSEXP, SEXP mu_inputSEXP, SEXP mu_dimSEXP, SEXP mu0_inputSEXP, SEXP mu0_dimSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type C(CSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_input(mu_inputSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mu_dim(mu_dimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0_input(mu0_inputSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mu0_dim(mu0_dimSEXP);
    Rcpp::traits::input_parameter< int >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(calibNoDist(Y, C, Z, mu_input, mu_dim, mu0_input, mu0_dim, ref));
    return rcpp_result_gen;
END_RCPP
}
// perturbedSNcpp
Rcpp::List perturbedSNcpp(arma::mat Y, arma::uvec C, Rcpp::List prior, Rcpp::List pmc, Rcpp::List state, Rcpp::List initParticles, bool init);
RcppExport SEXP _COMIX_perturbedSNcpp(SEXP YSEXP, SEXP CSEXP, SEXP priorSEXP, SEXP pmcSEXP, SEXP stateSEXP, SEXP initParticlesSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pmc(pmcSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type initParticles(initParticlesSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(perturbedSNcpp(Y, C, prior, pmc, state, initParticles, init));
    return rcpp_result_gen;
END_RCPP
}
// KL
double KL(arma::vec xi_1, arma::vec xi_2, arma::mat Omega_1, arma::mat Omega_2, arma::vec alpha_1, arma::vec alpha_2);
RcppExport SEXP _COMIX_KL(SEXP xi_1SEXP, SEXP xi_2SEXP, SEXP Omega_1SEXP, SEXP Omega_2SEXP, SEXP alpha_1SEXP, SEXP alpha_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type xi_1(xi_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xi_2(xi_2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega_1(Omega_1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega_2(Omega_2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha_1(alpha_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha_2(alpha_2SEXP);
    rcpp_result_gen = Rcpp::wrap(KL(xi_1, xi_2, Omega_1, Omega_2, alpha_1, alpha_2));
    return rcpp_result_gen;
END_RCPP
}
// relabel
Rcpp::List relabel(const Rcpp::List res);
RcppExport SEXP _COMIX_relabel(SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(relabel(res));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_COMIX_hungarian_cc", (DL_FUNC) &_COMIX_hungarian_cc, 1},
    {"_COMIX_calib", (DL_FUNC) &_COMIX_calib, 8},
    {"_COMIX_calibNoDist", (DL_FUNC) &_COMIX_calibNoDist, 8},
    {"_COMIX_perturbedSNcpp", (DL_FUNC) &_COMIX_perturbedSNcpp, 7},
    {"_COMIX_KL", (DL_FUNC) &_COMIX_KL, 6},
    {"_COMIX_relabel", (DL_FUNC) &_COMIX_relabel, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_COMIX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
