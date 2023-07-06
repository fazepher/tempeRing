#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


static double const logsqrt_tau = 0.5*log(arma::datum::tau);

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

// [[Rcpp::export]]
double lmvtnorm_temp_chol_cpp(const NumericVector& x, const NumericVector& mu,
                              double beta, const NumericMatrix& L_inv, double ldet_L_inv){
  int d = x.size();
  double cte = ldet_L_inv - d*(logsqrt_tau + 0.5*log(beta));
  return cte - beta*0.5*mahalanobis_chol_cpp(x, mu, L_inv);
}

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

// [[Rcpp::export]]
double dmvtnorm_temp_cpp(const NumericVector& x, const NumericVector& mu,
                         double beta = 1.0, double global_scale = 1.0,
                         const Nullable<NumericMatrix>& sigma = R_NilValue,
                         const Nullable<NumericMatrix>& L_inv = R_NilValue,
                         double ldet_L_inv = NA_REAL, bool log_p = false){

  double lp = lmvtnorm_temp_cpp(x, mu, beta, global_scale, sigma, L_inv, ldet_L_inv);
  return log_p ? lp : exp(lp);

}

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

// [[Rcpp::export]]
double lhatsn_cpp(const NumericVector& x, const NumericVector& m, double s, double a){

  NumericVector z = (x - m)/s;
  double cte = z.size()*(log(2) - log(s));
  return cte + sum(dnorm(z, 0.0, 1.0, true) + pnorm(a*z, 0.0, 1.0, true, true));

}

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

// [[Rcpp::export]]
double ulmixhatsn_temp_cpp(const NumericVector& x, double beta,
                           const NumericVector& w, const List& mu, const NumericVector& omega,
                           double alpha){

  return beta*lmixhatsn_cpp(x, w, mu, omega, alpha);

}

// [[Rcpp::export]]
List modAssignment_euclidean_cpp(const NumericVector& x, double beta,
                                 const NumericVector& w, const List& modes,
                                 const List& L_inv, const NumericVector& ldet_L_inv){

  int n_modes = w.size();
  NumericVector lp_m = log(w);
  for(int m = 0; m < n_modes; m++){
    lp_m[m] += lmvtnorm_temp_chol_cpp(x, modes[m], beta, L_inv[m], ldet_L_inv[m]);
  }
  int A = which_max(lp_m);

  return List::create(Named("A") = A + 1, Named("lP_j") = lp_m[A]);

}

// [[Rcpp::export]]
arma::vec lpsampler_cpp(const NumericVector& x_curr, double beta_max,
                        const NumericVector& w, const List& modes, const List& L){

  IntegerVector idx_modes = seq_along(w);
  int m = Rcpp::RcppArmadillo::sample(idx_modes, 1, false, w)[0] - 1;
  double temp_factor = exp(-0.5*log(beta_max));
  return rmvtnorm_chol_cpp(1, modes[m], L[m], temp_factor).col(0);

}

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
