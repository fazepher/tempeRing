#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


static double const logsqrt_tau = 0.5*log(arma::datum::tau);


////----Utilities----////



// [[Rcpp::export]]
double logsumexp_cpp(const NumericVector& lterms){

  int k = lterms.size();

  int imax = which_max(lterms);

  double acum = 0.0;
  for(int i = 0; i < k; i++){
    if(i == imax){
      acum += 1.0;
    }else{
      acum += exp(lterms[i] - lterms[imax]);
    }
  }

  return lterms[imax] + log(acum);

}

// [[Rcpp::export]]
double mahalanobis_chol_cpp(const NumericVector& x, const NumericVector& mu,
                            const NumericMatrix& L_inv, bool squared = true){

  int d = x.size();

  NumericVector z = x - mu;

  double maha_squared = 0.0;
  double li = 0.0;
  for(int i = 0; i < d; i++){

    // Iniciamos fila
    li = 0.0;
    // Triangular inferior, solo necesitamos multiplicar hasta la diagonal
    for(int j = 0; j <= i; j++){
      li += L_inv(i, j) * z[j];
    }
    // Vamos acumulando el producto punto
    maha_squared += li*li;

  }

  if(squared){
    return maha_squared;
  }
  return sqrt(maha_squared);

}



////----Distributions----////


//' @export
// [[Rcpp::export]]
double lmvtnorm_temp_chol_cpp(const NumericVector& x, const NumericVector& mu,
                              double beta, const NumericMatrix& L_inv, double ldet_L_inv){
  int d = x.size();
  double cte = ldet_L_inv - d*(logsqrt_tau + 0.5*log(beta));
  return cte - beta*0.5*mahalanobis_chol_cpp(x, mu, L_inv);
}
//' @export
// [[Rcpp::export]]
double lmvtnorm_temp_cpp(const NumericVector& x, const NumericVector& mu,
                         double beta = 1.0, double global_scale = 1.0,
                         const Nullable<NumericMatrix>& sigma = R_NilValue,
                         const Nullable<NumericMatrix>& L_inv = R_NilValue,
                         double ldet_L_inv = NA_REAL){

  // Obtenemos dimensión
  int d = x.size();

  if(L_inv.isNotNull()){ // Si el usuario da Cholesky, la usamos directamente

    NumericMatrix L_inv_(L_inv);
    double cte = ldet_L_inv - d*(logsqrt_tau + 0.5*log(beta));

    return cte - beta*0.5*mahalanobis_chol_cpp(x, mu, L_inv_);

  } else if(sigma.isNull()){ // Si no nos da ni Cholesky ni Sigma, es la identidad con escalamiento

    NumericVector z = x - mu;
    double scale_log_factor = 2*log(global_scale) - log(beta);
    double cte = d*(logsqrt_tau + scale_log_factor);

    return -0.5*(cte - exp(scale_log_factor)*sum(z*z));

  } else{ // Si nos da Sigma, calculamos Cholesky con arma; haciendo cálculo manual se evita copia a Rcpp

    NumericMatrix sigma_(sigma);
    arma::mat L_inv_ = inv(trimatl(chol(as<arma::mat>(sigma_), "lower")));

    ldet_L_inv = sum(log(L_inv_.diag()));
    double cte = ldet_L_inv - d*logsqrt_tau;

    arma::vec z = as<arma::vec>(x) - as<arma::vec>(mu);
    double maha_squared = 0.0;
    double li = 0.0;
    for(int i = 0; i < d; i++){
      li = 0.0;
      for(int j = 0; j <= i; j++){
        li += L_inv_(i, j) * z[j];
      }
      maha_squared += li*li;
    }

    return cte - 0.5*maha_squared;

  }

}
//' @export
// [[Rcpp::export]]
double dmvtnorm_temp_cpp(const NumericVector& x, const NumericVector& mu,
                         double beta = 1.0, double global_scale = 1.0,
                         const Nullable<NumericMatrix>& sigma = R_NilValue,
                         const Nullable<NumericMatrix>& L_inv = R_NilValue,
                         double ldet_L_inv = NA_REAL, bool log_p = false){

  double lp = lmvtnorm_temp_cpp(x, mu, beta, global_scale, sigma, L_inv, ldet_L_inv);
  return log_p ? lp : exp(lp);

}
//' @export
// [[Rcpp::export]]
arma::mat rmvtnorm_chol_cpp(int n, const arma::vec& mu, const arma::mat& L,
                            double scale_factor = 1.0){

  int d = mu.size();
  arma::mat x(d, n);
  arma::vec z(d);

  for(int k = 0; k < n; k++){

    x.col(k) += mu;
    z = scale_factor*rnorm(d); // Simulamos normales estándar con rescalamiento
    // Transformamos con matriz triangular inferior
    for(int i = 0; i < d; i++){
      for(int j = 0; j <= i; j++){
        x(i, k) += L(i, j) * z(j);
      }
    }

  }

  return x;

}

//' @export
// [[Rcpp::export]]
arma::mat rmvtnorm_temp_cpp(int n, const arma::vec& mu,
                            double beta = 1.0, double global_scale = 1.0,
                            const Nullable<NumericMatrix>& sigma = R_NilValue,
                            const Nullable<NumericMatrix>& L = R_NilValue){

  // Obtenemos dimensión y factor de escalamiento
  int d = mu.size();
  double scale_factor = beta == 1.0 ? global_scale : global_scale/sqrt(beta);

  if(L.isNotNull()){ // Si el usuario da Cholesky, la usamos directamente

    return rmvtnorm_chol_cpp(n, mu, as<arma::mat>(L), scale_factor);

  } else if(sigma.isNotNull()){ // Si nos dan Sigma, calculamos Cholesky con arma

    return rmvtnorm_chol_cpp(n, mu, chol(as<arma::mat>(sigma), "lower"), scale_factor);

  } else{ // Si no nos dan ni Cholesky ni Sigma, es la identidad con escalamiento

    arma::mat x(d, n);
    for(int i = 0; i < n; i++){
      x.col(i) += mu + global_scale*as<arma::colvec>(rnorm(d));
    }
    return x;

  }

}

//' @export
// [[Rcpp::export]]
double lhatsn_cpp(const NumericVector& x, const NumericVector& m, double s, double a){

  NumericVector z = (x - m)/s;
  double cte = z.size()*(log(2) - log(s));
  return cte + sum(dnorm(z, 0.0, 1.0, true) + pnorm(a*z, 0.0, 1.0, true, true));

}

//' @export
// [[Rcpp::export]]
double lmixhatsn_cpp(const NumericVector& x,
                     const NumericVector& w, const List& mu, const NumericVector& omega,
                     double alpha){

  int n_modes = w.size();
  NumericVector l_modes(n_modes);

  for(int m=0; m < n_modes; m++){
    l_modes[m] = log(w[m]) + lhatsn_cpp(x, mu[m], omega[m], alpha);
  }

  return logsumexp_cpp(l_modes);

}

//' @export
// [[Rcpp::export]]
double ulmixhatsn_temp_cpp(const NumericVector& x, double beta,
                           const NumericVector& w, const List& mu, const NumericVector& omega,
                           double alpha){

  return beta*lmixhatsn_cpp(x, w, mu, omega, alpha);

}

//' @export
// [[Rcpp::export]]
double ulmixhatsn_temp_cpp_alfas(const NumericVector& x, double beta,
                                 const NumericVector& w, const List& mu, const NumericVector& omega,
                                 const NumericVector& alpha){

  int n_modes = w.size();
  NumericVector l_modes(n_modes);

  for(int m=0; m < n_modes; m++){
    l_modes[m] = log(w[m]) + lhatsn_cpp(x, mu[m], omega[m], alpha[m]);
  }

  return beta*logsumexp_cpp(l_modes);

}



////----Metropolis-Hastings----////


//' @rdname mh_step
//'
// [[Rcpp::export]]
List mh_step_cpp(const NumericVector& x_curr, const NumericVector& x_prop,
                 double l_curr, double l_prop, double lq_c2p = 0.0, double lq_p2c = 0.0){

  int d = x_curr.size();
  NumericVector x_next(d);
  double l_next;

  double l_ratio = (l_prop - l_curr) + (lq_p2c - lq_c2p);
  double alpha = l_ratio > 0.0 ? 1.0 : exp(l_ratio);
  bool accepted = log(runif(1)[0]) <= l_ratio;

  if(accepted){
    x_next = x_prop;
    l_next = l_prop;
  }else{
    x_next = x_curr;
    l_next = l_curr;
  }


  return List::create(Named("x_next") = x_next,
                      Named("l_next") = l_next,
                      Named("accepted") = accepted,
                      Named("alpha") = alpha,
                      Named("l_ratio") = l_ratio);

}

//' @export
// [[Rcpp::export]]
List metropolis_step_cpp(const NumericVector& x_curr, const NumericVector& x_prop,
                         double l_curr, double l_prop){
  return mh_step_cpp(x_curr, x_prop, l_curr, l_prop);
}



////----ALPS----////



//' @export
// [[Rcpp::export]]
List modAssignment_cpp(const NumericVector& x, double beta,
                       const NumericVector& l_target_modes, const List& modes,
                       const List& L_inv, int n_modes){

  NumericVector half_maha(n_modes, 0.5);
  for(int m = 0; m < n_modes; m++){
    half_maha[m] *= mahalanobis_chol_cpp(x, modes[m], L_inv[m]);
  }
  int A_1 = which_max(l_target_modes - half_maha);
  NumericVector G_x_beta_m = l_target_modes - beta*half_maha;
  int A_beta = which_max(G_x_beta_m);

  return List::create(Named("A_beta") = A_beta + 1,
                      Named("A_1") = A_1 + 1,
                      Named("l_target_mu") = l_target_modes[A_beta],
                      Named("G_x_beta") = G_x_beta_m[A_beta]);

}

//' @export
// [[Rcpp::export]]
arma::vec lpsampler_cpp(const NumericVector& x_curr, double beta_max,
                        const NumericVector& w, const List& modes, const List& L){

  IntegerVector idx_modes = seq_along(w);
  int m = Rcpp::RcppArmadillo::sample(idx_modes, 1, false, w)[0] - 1;
  double temp_factor = exp(-0.5*log(beta_max));
  return rmvtnorm_chol_cpp(1, modes[m], L[m], temp_factor).col(0);

}

//' @export
// [[Rcpp::export]]
double lps_q_cpp(const NumericVector& x_curr, const NumericVector& x_prop,
                 double beta_max, const NumericVector& w, const List& modes,
                 const List& L_inv, const NumericVector& ldet_L_inv){

  NumericVector l_modes = log(w);
  for(int m = 0; m < w.size(); m++){
    l_modes[m] += lmvtnorm_temp_chol_cpp(x_prop, modes[m], beta_max, L_inv[m], ldet_L_inv[m]);
  }

  return logsumexp_cpp(l_modes);

}


////---- RJ-ALPS ----////

// [[Rcpp::export]]
List modAssignment_RJMCMC_cpp(const NumericVector& x, const NumericVector& x_tilde, double beta,
                              const NumericVector& l_target_modes,
                              const NumericVector& modes_idxs, const List& modes_tilde,
                              const List& L_inv, int n_modes){

  NumericVector half_maha(n_modes, arma::math::inf());
  for(int m = 0; m < n_modes; m++){
    if(modes_idxs[m] == x[0]){
      half_maha[m] = 0.5*mahalanobis_chol_cpp(x_tilde, modes_tilde[m], L_inv[m]);
    }
  }
  int A_1 = which_max(l_target_modes - half_maha);
  NumericVector G_x_beta_m = l_target_modes - beta*half_maha;
  int A_beta = which_max(G_x_beta_m);

  return List::create(Named("A_beta") = A_beta + 1,
                      Named("A_1") = A_1 + 1,
                      Named("l_target_mu") = l_target_modes[A_beta],
                      Named("G_x_beta") = G_x_beta_m[A_beta]);

}

// [[Rcpp::export]]
arma::mat rj_rmvtnorm_temp_cpp(int n, double r, const arma::vec& mu,
                               double beta = 1.0, double global_scale = 1.0,
                               const Nullable<NumericMatrix>& sigma = R_NilValue,
                               const Nullable<NumericMatrix>& L = R_NilValue){

  arma::mat x(mu.size(), n);
  x = rmvtnorm_temp_cpp(n, mu, beta, global_scale, sigma, L);
  arma::mat r_row(1, n);
  r_row.row(0) = r;
  arma::mat x_sampled(mu.size() + 1, n);
  x_sampled = join_cols(r_row, x);

  return x_sampled;

}

// [[Rcpp::export]]
arma::vec rj_lpsampler_cpp(const NumericVector& x_curr, double beta_max,
                           const NumericVector& w, const List& modes_tilde, const List& L){

  IntegerVector idx_modes = seq_along(w);

  int r_new = Rcpp::RcppArmadillo::sample(idx_modes, 1, false, w)[0];
  int m = r_new - 1;
  double temp_factor = exp(-0.5*log(beta_max));
  arma::vec param_mode = modes_tilde[m];
  arma::vec x_sample = rmvtnorm_chol_cpp(1, param_mode, L[m], temp_factor);

  arma::vec x_sample_full(x_sample.size() + 1);
  x_sample_full[0] = r_new;
  for(int i = 1; i < x_sample_full.size(); i++){
    x_sample_full[i] = x_sample[i-1];
  }

  return x_sample_full;

}

// [[Rcpp::export]]
double rj_lps_q_cpp(const NumericVector& x_curr, const NumericVector& x_prop,
                    double beta_max, const NumericVector& w, const List& modes,
                    const List& L_inv, const NumericVector& ldet_L_inv){

  int d_par = x_prop.size() - 1;
  IntegerVector idx_params = seq(1, d_par);
  LogicalVector valid_mode(w.size());
  NumericVector l_modes = log(w);
  for(int m = 0; m < w.size(); m++){
    NumericVector param_mode = modes[m];
    valid_mode[m] = x_prop[0] == param_mode[0];
    if(valid_mode[m]){
      l_modes[m] += lmvtnorm_temp_chol_cpp(x_prop[idx_params],
                                           param_mode[idx_params],
                                           beta_max, L_inv[m], ldet_L_inv[m]);
    }
  }

  return logsumexp_cpp(l_modes[valid_mode]);

}
