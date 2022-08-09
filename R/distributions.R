

# Basic Normal log densities and samplers
lnorm <- function(x, mean = 0, sd = 1){
  dnorm(x, mean, sd, log = TRUE)
}
rnorm <- stats::rnorm

lnorm_temp <- function(x, beta = 1, mean = 0, sd = 1){
  dnorm(x, mean, sd/sqrt(beta), log = TRUE)
}
rnorm_temp <- function(n, beta = 1, mean = 0, sd = 1){
  rnorm(n, mean, sd/sqrt(beta))
}

# Basic Multivariate normal log densities and samplers
lmvtnorm_temp <- function(x, beta = 1, mu = 0, sigma = NULL,
                          sigma_inv = NULL, logdet_sigma = NULL, d = NULL){

  # Dimension and number of points (determines vectorization behavior)
  if(is.matrix(x)){
    d <- d %||% ncol(x)
    n <- nrow(x)
  }else{
    d <- d %||% length(x)
    n <- 1
  }

  # If a single number provided for the mean, we expand to vector
  if(length(mu) == 1){
    mu <- rep(mu, d)
  }

  # If provided, we use the inverse; if not, by default sigma is the identity, otherwise we solve
  if(is.null(sigma_inv) && is.null(sigma)){
    sigma_inv <- diag(1, d)
    logdet_sigma <- 0
  }else{
    sigma_inv <- sigma_inv %||% solve(sigma)
    logdet_sigma <- logdet_sigma %||% as.vector(determinant(sigma)$modulus) #determinant returns log
  }

  # D2 = (x-mu)^T (sigma)^-1 (x-mu) is the squared Mahalanobis distance
  if(n == 1){
    D2 <- rowSums( (x - mu)%*%sigma_inv*(x - mu) )
  }else{
    D2 <- apply(x, 1, function(x) x - mu) |>
      apply(2, function(x_mu) rowSums( x_mu%*%sigma_inv*x_mu ))
  }

  if(beta == 1){
    return( -0.5*(d*log(2*pi) + logdet_sigma + D2) )
  }else{
    return( -0.5*(d*log(2*pi) - d*log(beta) + logdet_sigma + beta*D2) )
  }

}
rmvtnorm_temp <- function(n, beta = 1, mu, sigma = NULL,
                          LChol_sigma = NULL, d = NULL){

  d <- d %||% length(mu)
  LChol_sigma <- LChol_sigma %||% if(is.null(sigma)) diag(1,d) else t(chol(sigma))
  x <- replicate(n, as.vector( mu + (LChol_sigma/sqrt(beta)) %*% rnorm(d) ), simplify = TRUE)

  return(t(x))
}

lmvtnorm <- function(x, mu = 0, sigma = NULL, sigma_inv = NULL, logdet_sigma = NULL){
  lmvtnorm_temp(x = x, mu = mu, sigma = sigma, sigma_inv = sigma_inv, logdet_sigma = logdet_sigma)
}
rmvtnorm <- function(n, mu = 0, sigma = NULL, LChol_sigma = NULL){
  rmvtnorm_temp(n = n, mu = mu, sigma = sigma, LChol_sigma = LChol_sigma)
}

# General mixture distribution
lmix <- function(x, w, ldens, ..., shared_args = NULL){


  lk <- mapply(function(w, ...) log(w) + ldens(x, ...), w, ..., MoreArgs = shared_args)

  if(is.null(dim(lk))){
    return(matrixStats::logSumExp(lk))
  }else{
    return(apply(lk, 1, matrixStats::logSumExp))
  }


}
dmix <- function(x, w, ldens, ..., shared_args = NULL, log = FALSE){
  l <- lmix(x, w, ldens, ..., shared_args = shared_args)
  if(log){
    return(l)
  }
  return(exp(l))
}
rmix <- function(n, w, rdens, ..., shared_args = NULL){

  z <- sample.int(length(w), size = n, replace = TRUE, prob = w)
  x <- mapply(function(z,...) rdens(n = 1, ...), z, ... , MoreArgs = shared_args)

  return(x)
}

# General tempered unnormalized densities
ul_temp <- function(x, beta = 1, ldens, ...){
  beta*ldens(x, ...)
}
ulmix_temp <- function(x, beta = 1, w, ldens, ..., shared_args = NULL){
  # debug_aux <- match.call()
  # print(debug_aux)
  beta*lmix(x, w, ldens, ..., shared_args = shared_args)
}
lmix_temp <- function(x, beta = 1, w, ldens, ..., shared_args = NULL, z = NULL){

  aux_fun <- function(x, aux_beta = beta, aux_w = w,  aux_ldens = ldens, ...,
                      aux_shared_args = shared_args, z = z){

    ul <- ulmix_temp(x, aux_beta, aux_w, aux_ldens, ..., shared_args = aux_shared_args)

    if(is.null(z)){
      return(exp(ul))
    }else{
      return(ul)
    }

  }

  if(is.null(z)){
    z <- integrate(aux_fun, aux_beta = beta, aux_w = w,
                   aux_ldens = ldens, ..., aux_shared_args = shared_args, z = z,
                   lower = -Inf, upper = Inf)$value
  }

  return(aux_fun(x, ..., z = z) - log(z))

}

# Normal Mixture
lmix_norm <- function(x, w, mean, sd = NULL, shared_sd = 1){
  if(is.null(sd)){
    lmix(x, w, lnorm, mean, shared_args = list(sd = shared_sd))
  }else{
    lmix(x, w, lnorm, mean, sd)
  }
}
dmix_norm <- function(x, w, mean, sd = NULL, shared_sd = 1, log = FALSE){
  if(is.null(sd)){
    dmix(x, w, lnorm, mean, shared_args = list(sd = shared_sd), log = log)
  }else{
    dmix(x, w, lnorm, mean, sd, log = log)
  }
}
rmix_norm <- function(n, w, mean, sd = NULL, shared_sd = 1){
  if(is.null(sd)){
    rmix(n, w, rnorm, mean, shared_args = list(sd = shared_sd))
  }else{
    rmix(n, w, rnorm, mean, sd)
  }

}

ulmix_norm_temp <- function(x, beta = 1, w, mean, sd = NULL, shared_sd = 1){
  ul_temp(x, beta, lmix_norm, w, mean, sd, shared_sd)
}
# lmix_norm_temp <- function(x, beta = 1, w, mean, sd, shared_sd){
#   z, integrate_args = list(lower = -Inf, upper =))
# }

ulNorm_mix_temp <- function(x, beta = 1, w = c(0.5, 0.5), mu = c(-5, 5), sd = c(1, 1)){
  ulDens_mix_temp(x, w, ldens = lnorm_temp, ..., common_args = list(beta = beta))
}
nlmix_norm_temp <- function(x, beta = 1, w = c(0.5, 0.5), ldens, ..., common_args = NULL,
                            z = NULL, integrate_args = list(lower = -Inf, upper = Inf)){

  lbeta <- function(y){ ulmix_dens(y, w, ldens, ..., common_args = list(beta = beta)) }
  z <- z %||% do.call(integrate, c(list(f = lbeta), integrate_args))$value
  return(lbeta(x)/z)

}
#ulmixmvtnorm_temp
