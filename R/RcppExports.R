# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

logsumexp_cpp <- function(lterms) {
    .Call(`_tempeRing_logsumexp_cpp`, lterms)
}

mahalanobis_chol_cpp <- function(x, mu, L_inv, squared = TRUE) {
    .Call(`_tempeRing_mahalanobis_chol_cpp`, x, mu, L_inv, squared)
}

#' @export
lmvtnorm_temp_chol_cpp <- function(x, mu, beta, L_inv, ldet_L_inv) {
    .Call(`_tempeRing_lmvtnorm_temp_chol_cpp`, x, mu, beta, L_inv, ldet_L_inv)
}

#' @export
lmvtnorm_temp_cpp <- function(x, mu, beta = 1.0, global_scale = 1.0, sigma = NULL, L_inv = NULL, ldet_L_inv = NA_real_) {
    .Call(`_tempeRing_lmvtnorm_temp_cpp`, x, mu, beta, global_scale, sigma, L_inv, ldet_L_inv)
}

#' @export
dmvtnorm_temp_cpp <- function(x, mu, beta = 1.0, global_scale = 1.0, sigma = NULL, L_inv = NULL, ldet_L_inv = NA_real_, log_p = FALSE) {
    .Call(`_tempeRing_dmvtnorm_temp_cpp`, x, mu, beta, global_scale, sigma, L_inv, ldet_L_inv, log_p)
}

#' @export
rmvtnorm_chol_cpp <- function(n, mu, L, scale_factor = 1.0) {
    .Call(`_tempeRing_rmvtnorm_chol_cpp`, n, mu, L, scale_factor)
}

#' @export
rmvtnorm_temp_cpp <- function(n, mu, beta = 1.0, global_scale = 1.0, sigma = NULL, L = NULL) {
    .Call(`_tempeRing_rmvtnorm_temp_cpp`, n, mu, beta, global_scale, sigma, L)
}

#' @export
lhatsn_cpp <- function(x, m, s, a) {
    .Call(`_tempeRing_lhatsn_cpp`, x, m, s, a)
}

#' @export
lmixhatsn_cpp <- function(x, w, mu, omega, alpha) {
    .Call(`_tempeRing_lmixhatsn_cpp`, x, w, mu, omega, alpha)
}

#' @export
ulmixhatsn_temp_cpp <- function(x, beta, w, mu, omega, alpha) {
    .Call(`_tempeRing_ulmixhatsn_temp_cpp`, x, beta, w, mu, omega, alpha)
}

#' @export
ulmixhatsn_temp_cpp_alfas <- function(x, beta, w, mu, omega, alpha) {
    .Call(`_tempeRing_ulmixhatsn_temp_cpp_alfas`, x, beta, w, mu, omega, alpha)
}

#' @rdname mh_step
#'
mh_step_cpp <- function(x_curr, x_prop, l_curr, l_prop, lq_c2p = 0.0, lq_p2c = 0.0) {
    .Call(`_tempeRing_mh_step_cpp`, x_curr, x_prop, l_curr, l_prop, lq_c2p, lq_p2c)
}

#' @export
metropolis_step_cpp <- function(x_curr, x_prop, l_curr, l_prop) {
    .Call(`_tempeRing_metropolis_step_cpp`, x_curr, x_prop, l_curr, l_prop)
}

#' @export
modAssignment_cpp <- function(x, beta, l_target_modes, modes, L_inv, n_modes) {
    .Call(`_tempeRing_modAssignment_cpp`, x, beta, l_target_modes, modes, L_inv, n_modes)
}

#' @export
lpsampler_cpp <- function(x_curr, beta_max, w, modes, L) {
    .Call(`_tempeRing_lpsampler_cpp`, x_curr, beta_max, w, modes, L)
}

#' @export
lps_q_cpp <- function(x_curr, x_prop, beta_max, w, modes, L_inv, ldet_L_inv) {
    .Call(`_tempeRing_lps_q_cpp`, x_curr, x_prop, beta_max, w, modes, L_inv, ldet_L_inv)
}

modAssignment_RJMCMC_cpp <- function(x, x_tilde, beta, l_target_modes, modes_idxs, modes_tilde, L_inv, n_modes) {
    .Call(`_tempeRing_modAssignment_RJMCMC_cpp`, x, x_tilde, beta, l_target_modes, modes_idxs, modes_tilde, L_inv, n_modes)
}

rj_rmvtnorm_temp_cpp <- function(n, r, mu, beta = 1.0, global_scale = 1.0, sigma = NULL, L = NULL) {
    .Call(`_tempeRing_rj_rmvtnorm_temp_cpp`, n, r, mu, beta, global_scale, sigma, L)
}

rj_lpsampler_cpp <- function(x_curr, beta_max, w, modes_tilde, L) {
    .Call(`_tempeRing_rj_lpsampler_cpp`, x_curr, beta_max, w, modes_tilde, L)
}

rj_lps_q_cpp <- function(x_curr, x_prop, beta_max, w, modes, L_inv, ldet_L_inv) {
    .Call(`_tempeRing_rj_lps_q_cpp`, x_curr, x_prop, beta_max, w, modes, L_inv, ldet_L_inv)
}

