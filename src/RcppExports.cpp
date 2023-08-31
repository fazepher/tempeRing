// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logsumexp_cpp
double logsumexp_cpp(const NumericVector& lterms);
RcppExport SEXP _tempeRing_logsumexp_cpp(SEXP ltermsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type lterms(ltermsSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp_cpp(lterms));
    return rcpp_result_gen;
END_RCPP
}
// mahalanobis_chol_cpp
double mahalanobis_chol_cpp(const NumericVector& x, const NumericVector& mu, const NumericMatrix& L_inv, bool squared);
RcppExport SEXP _tempeRing_mahalanobis_chol_cpp(SEXP xSEXP, SEXP muSEXP, SEXP L_invSEXP, SEXP squaredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type L_inv(L_invSEXP);
    Rcpp::traits::input_parameter< bool >::type squared(squaredSEXP);
    rcpp_result_gen = Rcpp::wrap(mahalanobis_chol_cpp(x, mu, L_inv, squared));
    return rcpp_result_gen;
END_RCPP
}
// lmvtnorm_temp_chol_cpp
double lmvtnorm_temp_chol_cpp(const NumericVector& x, const NumericVector& mu, double beta, const NumericMatrix& L_inv, double ldet_L_inv);
RcppExport SEXP _tempeRing_lmvtnorm_temp_chol_cpp(SEXP xSEXP, SEXP muSEXP, SEXP betaSEXP, SEXP L_invSEXP, SEXP ldet_L_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type L_inv(L_invSEXP);
    Rcpp::traits::input_parameter< double >::type ldet_L_inv(ldet_L_invSEXP);
    rcpp_result_gen = Rcpp::wrap(lmvtnorm_temp_chol_cpp(x, mu, beta, L_inv, ldet_L_inv));
    return rcpp_result_gen;
END_RCPP
}
// lmvtnorm_temp_cpp
double lmvtnorm_temp_cpp(const NumericVector& x, const NumericVector& mu, double beta, double global_scale, const Nullable<NumericMatrix>& sigma, const Nullable<NumericMatrix>& L_inv, double ldet_L_inv);
RcppExport SEXP _tempeRing_lmvtnorm_temp_cpp(SEXP xSEXP, SEXP muSEXP, SEXP betaSEXP, SEXP global_scaleSEXP, SEXP sigmaSEXP, SEXP L_invSEXP, SEXP ldet_L_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type global_scale(global_scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type L_inv(L_invSEXP);
    Rcpp::traits::input_parameter< double >::type ldet_L_inv(ldet_L_invSEXP);
    rcpp_result_gen = Rcpp::wrap(lmvtnorm_temp_cpp(x, mu, beta, global_scale, sigma, L_inv, ldet_L_inv));
    return rcpp_result_gen;
END_RCPP
}
// dmvtnorm_temp_cpp
double dmvtnorm_temp_cpp(const NumericVector& x, const NumericVector& mu, double beta, double global_scale, const Nullable<NumericMatrix>& sigma, const Nullable<NumericMatrix>& L_inv, double ldet_L_inv, bool log_p);
RcppExport SEXP _tempeRing_dmvtnorm_temp_cpp(SEXP xSEXP, SEXP muSEXP, SEXP betaSEXP, SEXP global_scaleSEXP, SEXP sigmaSEXP, SEXP L_invSEXP, SEXP ldet_L_invSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type global_scale(global_scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type L_inv(L_invSEXP);
    Rcpp::traits::input_parameter< double >::type ldet_L_inv(ldet_L_invSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvtnorm_temp_cpp(x, mu, beta, global_scale, sigma, L_inv, ldet_L_inv, log_p));
    return rcpp_result_gen;
END_RCPP
}
// rmvtnorm_chol_cpp
arma::mat rmvtnorm_chol_cpp(int n, const arma::vec& mu, const arma::mat& L, double scale_factor);
RcppExport SEXP _tempeRing_rmvtnorm_chol_cpp(SEXP nSEXP, SEXP muSEXP, SEXP LSEXP, SEXP scale_factorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type scale_factor(scale_factorSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvtnorm_chol_cpp(n, mu, L, scale_factor));
    return rcpp_result_gen;
END_RCPP
}
// rmvtnorm_temp_cpp
arma::mat rmvtnorm_temp_cpp(int n, const arma::vec& mu, double beta, double global_scale, const Nullable<NumericMatrix>& sigma, const Nullable<NumericMatrix>& L);
RcppExport SEXP _tempeRing_rmvtnorm_temp_cpp(SEXP nSEXP, SEXP muSEXP, SEXP betaSEXP, SEXP global_scaleSEXP, SEXP sigmaSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type global_scale(global_scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvtnorm_temp_cpp(n, mu, beta, global_scale, sigma, L));
    return rcpp_result_gen;
END_RCPP
}
// lhatsn_cpp
double lhatsn_cpp(const NumericVector& x, const NumericVector& m, double s, double a);
RcppExport SEXP _tempeRing_lhatsn_cpp(SEXP xSEXP, SEXP mSEXP, SEXP sSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(lhatsn_cpp(x, m, s, a));
    return rcpp_result_gen;
END_RCPP
}
// lmixhatsn_cpp
double lmixhatsn_cpp(const NumericVector& x, const NumericVector& w, const List& mu, const NumericVector& omega, double alpha);
RcppExport SEXP _tempeRing_lmixhatsn_cpp(SEXP xSEXP, SEXP wSEXP, SEXP muSEXP, SEXP omegaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const List& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(lmixhatsn_cpp(x, w, mu, omega, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ulmixhatsn_temp_cpp
double ulmixhatsn_temp_cpp(const NumericVector& x, double beta, const NumericVector& w, const List& mu, const NumericVector& omega, double alpha);
RcppExport SEXP _tempeRing_ulmixhatsn_temp_cpp(SEXP xSEXP, SEXP betaSEXP, SEXP wSEXP, SEXP muSEXP, SEXP omegaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const List& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ulmixhatsn_temp_cpp(x, beta, w, mu, omega, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ulmixhatsn_temp_cpp_alfas
double ulmixhatsn_temp_cpp_alfas(const NumericVector& x, double beta, const NumericVector& w, const List& mu, const NumericVector& omega, const NumericVector& alpha);
RcppExport SEXP _tempeRing_ulmixhatsn_temp_cpp_alfas(SEXP xSEXP, SEXP betaSEXP, SEXP wSEXP, SEXP muSEXP, SEXP omegaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const List& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ulmixhatsn_temp_cpp_alfas(x, beta, w, mu, omega, alpha));
    return rcpp_result_gen;
END_RCPP
}
// mh_step_cpp
List mh_step_cpp(const NumericVector& x_curr, const NumericVector& x_prop, double l_curr, double l_prop, double lq_c2p, double lq_p2c);
RcppExport SEXP _tempeRing_mh_step_cpp(SEXP x_currSEXP, SEXP x_propSEXP, SEXP l_currSEXP, SEXP l_propSEXP, SEXP lq_c2pSEXP, SEXP lq_p2cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_curr(x_currSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_prop(x_propSEXP);
    Rcpp::traits::input_parameter< double >::type l_curr(l_currSEXP);
    Rcpp::traits::input_parameter< double >::type l_prop(l_propSEXP);
    Rcpp::traits::input_parameter< double >::type lq_c2p(lq_c2pSEXP);
    Rcpp::traits::input_parameter< double >::type lq_p2c(lq_p2cSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_step_cpp(x_curr, x_prop, l_curr, l_prop, lq_c2p, lq_p2c));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_step_cpp
List metropolis_step_cpp(const NumericVector& x_curr, const NumericVector& x_prop, double l_curr, double l_prop);
RcppExport SEXP _tempeRing_metropolis_step_cpp(SEXP x_currSEXP, SEXP x_propSEXP, SEXP l_currSEXP, SEXP l_propSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_curr(x_currSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_prop(x_propSEXP);
    Rcpp::traits::input_parameter< double >::type l_curr(l_currSEXP);
    Rcpp::traits::input_parameter< double >::type l_prop(l_propSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_step_cpp(x_curr, x_prop, l_curr, l_prop));
    return rcpp_result_gen;
END_RCPP
}
// modAssignment_cpp
List modAssignment_cpp(const NumericVector& x, double beta, const NumericVector& l_target_modes, const List& modes, const List& L_inv, int n_modes);
RcppExport SEXP _tempeRing_modAssignment_cpp(SEXP xSEXP, SEXP betaSEXP, SEXP l_target_modesSEXP, SEXP modesSEXP, SEXP L_invSEXP, SEXP n_modesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type l_target_modes(l_target_modesSEXP);
    Rcpp::traits::input_parameter< const List& >::type modes(modesSEXP);
    Rcpp::traits::input_parameter< const List& >::type L_inv(L_invSEXP);
    Rcpp::traits::input_parameter< int >::type n_modes(n_modesSEXP);
    rcpp_result_gen = Rcpp::wrap(modAssignment_cpp(x, beta, l_target_modes, modes, L_inv, n_modes));
    return rcpp_result_gen;
END_RCPP
}
// lpsampler_cpp
arma::vec lpsampler_cpp(const NumericVector& x_curr, double beta_max, const NumericVector& w, const List& modes, const List& L);
RcppExport SEXP _tempeRing_lpsampler_cpp(SEXP x_currSEXP, SEXP beta_maxSEXP, SEXP wSEXP, SEXP modesSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_curr(x_currSEXP);
    Rcpp::traits::input_parameter< double >::type beta_max(beta_maxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const List& >::type modes(modesSEXP);
    Rcpp::traits::input_parameter< const List& >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(lpsampler_cpp(x_curr, beta_max, w, modes, L));
    return rcpp_result_gen;
END_RCPP
}
// lps_q_cpp
double lps_q_cpp(const NumericVector& x_curr, const NumericVector& x_prop, double beta_max, const NumericVector& w, const List& modes, const List& L_inv, const NumericVector& ldet_L_inv);
RcppExport SEXP _tempeRing_lps_q_cpp(SEXP x_currSEXP, SEXP x_propSEXP, SEXP beta_maxSEXP, SEXP wSEXP, SEXP modesSEXP, SEXP L_invSEXP, SEXP ldet_L_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x_curr(x_currSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_prop(x_propSEXP);
    Rcpp::traits::input_parameter< double >::type beta_max(beta_maxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const List& >::type modes(modesSEXP);
    Rcpp::traits::input_parameter< const List& >::type L_inv(L_invSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type ldet_L_inv(ldet_L_invSEXP);
    rcpp_result_gen = Rcpp::wrap(lps_q_cpp(x_curr, x_prop, beta_max, w, modes, L_inv, ldet_L_inv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tempeRing_logsumexp_cpp", (DL_FUNC) &_tempeRing_logsumexp_cpp, 1},
    {"_tempeRing_mahalanobis_chol_cpp", (DL_FUNC) &_tempeRing_mahalanobis_chol_cpp, 4},
    {"_tempeRing_lmvtnorm_temp_chol_cpp", (DL_FUNC) &_tempeRing_lmvtnorm_temp_chol_cpp, 5},
    {"_tempeRing_lmvtnorm_temp_cpp", (DL_FUNC) &_tempeRing_lmvtnorm_temp_cpp, 7},
    {"_tempeRing_dmvtnorm_temp_cpp", (DL_FUNC) &_tempeRing_dmvtnorm_temp_cpp, 8},
    {"_tempeRing_rmvtnorm_chol_cpp", (DL_FUNC) &_tempeRing_rmvtnorm_chol_cpp, 4},
    {"_tempeRing_rmvtnorm_temp_cpp", (DL_FUNC) &_tempeRing_rmvtnorm_temp_cpp, 6},
    {"_tempeRing_lhatsn_cpp", (DL_FUNC) &_tempeRing_lhatsn_cpp, 4},
    {"_tempeRing_lmixhatsn_cpp", (DL_FUNC) &_tempeRing_lmixhatsn_cpp, 5},
    {"_tempeRing_ulmixhatsn_temp_cpp", (DL_FUNC) &_tempeRing_ulmixhatsn_temp_cpp, 6},
    {"_tempeRing_ulmixhatsn_temp_cpp_alfas", (DL_FUNC) &_tempeRing_ulmixhatsn_temp_cpp_alfas, 6},
    {"_tempeRing_mh_step_cpp", (DL_FUNC) &_tempeRing_mh_step_cpp, 6},
    {"_tempeRing_metropolis_step_cpp", (DL_FUNC) &_tempeRing_metropolis_step_cpp, 4},
    {"_tempeRing_modAssignment_cpp", (DL_FUNC) &_tempeRing_modAssignment_cpp, 6},
    {"_tempeRing_lpsampler_cpp", (DL_FUNC) &_tempeRing_lpsampler_cpp, 5},
    {"_tempeRing_lps_q_cpp", (DL_FUNC) &_tempeRing_lps_q_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_tempeRing(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
