
##### Basic Normal log densities and samplers #####

#' Tempered Normal log densities and sampler
#'
#' @description
#' These are basic wrappers around `stats::dnorm()` and `stats::rnorm()` to provide shortcut evaluation
#' of the log-density, density and random generation of a univariate tempered normal
#' with inverse temperature equal to `beta`, mean equal to `mean` and **standard deviation** equal to `sd`.
#'
#' When `beta` = 1, we recover a regular normal density.
#'
#' @details
#' Tempering a distribution means raising its density to a power \eqn{\beta>0},
#' known as inverse temperature. Equivalently, we multiply the log-density by \eqn{\beta}:
#' \deqn{f_\beta(x) = f(x)^\beta}
#' \deqn{l_\beta(x) = \beta l(x)}
#' Consider a univariate normal random variable centered at \eqn{\mu} with
#' standard deviation \eqn{\sigma}, where \eqn{cte} represents the normalizing constant
#' \deqn{X ~ N(\mu, \sigma)}
#' \deqn{l(x) = -(x-\mu)^2 / 2\sigma^2 + cte}
#' Its tempered version is equivalent to rescaling with new standard deviation \eqn{\sigma/\sqrt\beta}
#' and keeping the same mean parameter:
#' \deqn{l_\beta(x) = \beta l(x) = -\beta(x-\mu)^2 / 2\sigma^2 + cte'}
#' \deqn{X|\beta ~ N(\mu, \sigma/\sqrt\beta)}
#'
#'
#' @param x vector of quantiles.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param beta inverse temperature parameter \eqn{\beta > 0}.
#' @param n number of observations.
#'
#' @returns
#' The preffix `l` stands for log-density, `d` for density, and `r` for sampling.
#' `lnorm` gives the log-density of a regular normal without tempering (i.e. \eqn{\beta = 1}).
#' See [stats::dnorm()] for more information.
#'
#' @seealso [lmvtnorm_temp()], [stats::dnorm()]
#'
#' @export
#'
#' @examples
#'
#' lnorm_temp(x = 0, beta = 0.5, sd = 1)
#' lnorm(0, sd = 1/sqrt(0.5))
#'
#' lnorm(x = 0)
#' dnorm(x = 0, log = TRUE)
#'
#' dnorm_temp(x = 5, beta = 1)
#' dnorm(x = 5)
#'
#' rnorm_temp(n = 1000, mean = 100, sd = 1) |> hist()
#' rnorm_temp(n = 1000, beta = 0.1, mean = 100, sd = 1) |> hist()
#'
#' # The functions inherit vectorization,
#' # so can be used for example with ggplot2::geom_function()
#' # to show the flattening effect of tempering
#'
#'if(requireNamespace("ggplot2", quietly = TRUE)){
#' ggplot() +
#'   geom_function(fun = dnorm,
#'                 color = "gray65") +
#'   geom_function(fun = dnorm_temp,
#'                 args = list(beta = 0.75),
#'                 color = "steelblue4") +
#'   geom_function(fun = dnorm_temp,
#'                 args = list(beta = 0.5),
#'                 color = "darkcyan") +
#'   geom_function(fun = dnorm_temp,
#'                 args = list(beta = 0.25),
#'                 color = "blueviolet") +
#'   geom_function(fun = dnorm_temp,
#'                 args = list(beta = 0.25^2),
#'                 color = "maroon4") +
#'   xlim(-6,6) +
#'   theme_classic()
#'}
#'
#'
lnorm_temp <- function(x, beta = 1, mean = 0, sd = 1){
  dnorm(x, mean, sd/sqrt(beta), log = TRUE)
}

#' @rdname lnorm_temp
#'
#' @param log For `dnorm_temp`, whether to return the log-density or not (the default).
#'
#' @export
#'
dnorm_temp <- function(x, beta = 1, mean = 0, sd = 1, log = FALSE){
  dnorm(x, mean, sd/sqrt(beta), log)
}

#' @rdname lnorm_temp
#'
#' @export
#'
rnorm_temp <- function(n, beta = 1, mean = 0, sd = 1){
  rnorm(n, mean, sd/sqrt(beta))
}

#' @rdname lnorm_temp
#'
#' @export
#'
lnorm <- function(x, mean = 0, sd = 1){
  dnorm(x, mean, sd, log = TRUE)
}


##### Basic Multivariate normal log densities and samplers ####

#' Multivariate Tempered Normal log densities and sampler
#'
#' @description
#'
#' Log-density, density and random generation of a multivariate tempered normal
#' with inverse temperature equal to `beta`, mean equal to `mu` and **covariance** matrix
#' equal to `sigma`.
#'
#' `lmvtnorm`, `dmvtnorm` and `rmvtnorm` are equivalent to the tempered versions with `beta`=1.
#'
#' @details
#' Tempering a distribution means raising its density to a power \eqn{\beta>0},
#' known as inverse temperature. Equivalently, we multiply the log-density by \eqn{\beta}:
#' \deqn{f_\beta(x) = f(x)^\beta}
#' \deqn{l_\beta(x) = \beta l(x)}
#' Consider a \eqn{d}-dimensional multivariate normal random variable centered at \eqn{\mu} with
#' covariance matrix \eqn{\Sigma}, where \eqn{cte} represents the the normalizing constant
#' \deqn{X ~ MVN(\mu, \Sigma)}
#' \deqn{l(x) = -0.5(x-\mu)^T\Sigma^-1(x-\mu) + cte}
#' Its tempered version is equivalent to rescaling with new covariance matrix \eqn{\Sigma/\beta}
#' and keeping the same mean parameter:
#' \deqn{l_\beta(x) = -\beta 0.5 (x-\mu)^T \Sigma^-1 (x-\mu)  + cte' }
#' \deqn{l_\beta(x) = -0.5 (x-\mu)^T (\Sigma/\beta)^-1 (x-\mu) + cte'}
#' \deqn{X|\beta ~ MVN(\mu, \Sigma/\beta)}
#'
#' Now, the multivariate normal density only depends on \eqn{\Sigma} through its inverse
#' \eqn{\Sigma^-1} in the kernel term and determinant in the normalizing constant
#' \deqn{cte = -0.5(d log(2\pi) + log( det(\Sigma) )}
#' For this reason, instead of providing the covariance matrix `sigma`, the user can provide both via
#' `sigma_inv` and `logdet_sigma`, saving the functions the need to compute them under the hood.
#'
#' Another way of thinking of a Multivariate Normal random variable is to consider a Cholesky
#' decomposition of \eqn{\Sigma = LL^T}, whose lower triangular component \eqn{L} acts as
#' a linear transformation of \eqn{d} independent univariate standard normal variables in a vector
#' \deqn{X = \mu + L [z_1, ..., z_d]^T}
#' This is the way `rmvtnorm_temp` and `rmvtnorm` generate samples, via `stats::rnorm`.
#' Thus, if known, the user may provide `LChol_sigma` instead of `sigma`.
#'
#'
#' @param x A quantile vector or a matrix of quantile vectors (by rows).
#' @param beta Inverse temperature parameter \eqn{\beta > 0}.
#' @param mu Mean vector. If a single value is provided, it is expanded via `rep(mu,d)`.
#' @param sigma Covariance matrix, by default it is taken to be the identity.
#' @param sigma_inv (Optional) Inverse of the covariance matrix, see Details.
#' @param logdet_sigma (Optional) Logarithm of the determinant of the covariance matrix, see Details.
#' @param d Dimension of `x`, if NULL (the default) it is taken to be `ncol(x)` or `length(x)`
#' as appropriate;  for random generation, it would be taken as `length(mu)`.
#'
#' @returns
#' The preffix `l` stands for log-density, `d` for density, and `r` for sampling.
#' `lnorm` gives the log-density of a regular normal without tempering (i.e. \eqn{\beta = 1}).
#' See [stats::dnorm()] for more information. When sampling `n`> 1 realizations,
#' the resulting matrix has `n` rows and `d` columns.
#'
#' @examples
#'
#' Sigma <- matrix(c(1,0.5,0.5,1),2)
#' Sigma_inv <- solve(Sigma)
#' LogDet_Sigma <- log(det(Sigma))
#'
#' lmvtnorm_temp(x = c(0,0), beta = 0.5, sigma = Sigma)
#' lmvtnorm_temp(x = c(0,0), beta = 0.5,
#'               sigma_inv = Sigma_inv, logdet_sigma = LogDet_Sigma)
#' lmvtnorm(x = c(0,0), sigma = Sigma/0.5)
#' dmvtnorm(x = c(0,0), sigma = Sigma/0.5) |> log()
#'
#' L <- t(chol(Sigma))
#' rmvtnorm_temp(n = 1000, sigma = Sigma) |> plot()
#' rmvtnorm_temp(n = 1000, LChol_sigma = L) |> plot()
#'
#' rmvtnorm_temp(n = 1, mu = rep(5,3))
#' rmvtnorm_temp(n = 2, mu = rep(5,3))
#' rmvtnorm_temp(n = 2, mu = 5, d = 3)
#'
#' @export
#'
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
#' @rdname lmvtnorm_temp
#'
#' @param log For `dmvtnorm_temp` and `dmvtnorm`, whether to return the log-density or not (the default).
#'
#' @export
#'
dmvtnorm_temp <- function(x, beta = 1, mu = 0, sigma = NULL,
                          sigma_inv = NULL, logdet_sigma = NULL, d = NULL, log = FALSE){
  l <- lmvtnorm_temp(x, beta, mu, sigma, sigma_inv, logdet_sigma, d)
  if(log) return(l) else return(exp(l))
}
#' @rdname lmvtnorm_temp
#'
#' @param LChol_sigma (Optional) Lower triangular component of the Cholesky decomposition of `sigma`,
#' see Details.
#'
#' @export
#'
rmvtnorm_temp <- function(n, beta = 1, mu = rep(0,2), sigma = NULL,
                          LChol_sigma = NULL, d = NULL){

  d <- d %||% length(mu)
  LChol_sigma <- LChol_sigma %||% if(is.null(sigma)) diag(1,d) else t(chol(sigma))
  x <- replicate(n, as.vector( mu + (LChol_sigma/sqrt(beta)) %*% rnorm(d) ), simplify = TRUE)

  return(t(x))
}

#' @rdname lmvtnorm_temp
#'
#' @export
#'
lmvtnorm <- function(x, mu = 0, sigma = NULL, sigma_inv = NULL, logdet_sigma = NULL){
  lmvtnorm_temp(x = x, mu = mu, sigma = sigma, sigma_inv = sigma_inv, logdet_sigma = logdet_sigma)
}
#' @rdname lmvtnorm_temp
#'
#' @export
#'
dmvtnorm <- function(x, mu = 0, sigma = NULL, sigma_inv = NULL, logdet_sigma = NULL, log = FALSE){
  dmvtnorm_temp(x = x, mu = mu, sigma = sigma,
                sigma_inv = sigma_inv, logdet_sigma = logdet_sigma,
                log = log)
}
#' @rdname lmvtnorm_temp
#'
#' @export
#'
rmvtnorm <- function(n, mu = 0, sigma = NULL, LChol_sigma = NULL){
  rmvtnorm_temp(n = n, mu = mu, sigma = sigma, LChol_sigma = LChol_sigma)
}

####  General mixture distribution ####

#' General Finite Mixture Distribution
#'
#' @export
#'
lmix <- function(x, w, ldens, ..., shared_args = NULL){

  lk <- mapply(function(w, ...) log(w) + ldens(x, ...), w, ..., MoreArgs = shared_args)

  if(is.null(dim(lk))){
    return(logSumExp(lk))
  }else{
    return(apply(lk, 1, logSumExp))
  }

}
#' @rdname lmix
#'
#' @export
#'
dmix <- function(x, w, ldens, ..., shared_args = NULL, log = FALSE){
  l <- lmix(x, w, ldens, ..., shared_args = shared_args)
  if(log) return(l) else return(exp(l))
}
#' @rdname lmix
#'
#' @export
#'
rmix <- function(n, w, rdens, ..., shared_args = NULL){
  z <- sample.int(length(w), size = n, replace = TRUE, prob = w)
  x <- mapply(function(z,...) rdens(n = 1, ...), z, ... , MoreArgs = shared_args)
  return(x)
}

#### General tempered densities ####

#' General Tempered density generator
#'
#' @export
#'
ul_temp <- function(x, beta = 1, ldens, ...) beta*ldens(x, ...)


#' General Tempered Mixture generator
#'
#' @export
#'
ulmix_temp <- function(x, beta = 1, w, ldens, ..., shared_args = NULL){
  beta*lmix(x, w, ldens, ..., shared_args = shared_args)
}

lmix_temp <- function(x, beta = 1, w, ldens, ..., shared_args = NULL, z = NULL){

  aux_fun <- function(x, aux_beta = beta, aux_w = w,  aux_ldens = ldens, ...,
                      aux_shared_args = shared_args, z = z){

    ul <- ulmix_temp(x, aux_beta, aux_w, aux_ldens, ..., shared_args = aux_shared_args)

    if(is.null(z)) return(exp(ul)) else return(ul)

  }

  if(is.null(z)){
    z <- stats::integrate(aux_fun, aux_beta = beta, aux_w = w,
                          aux_ldens = ldens, ..., aux_shared_args = shared_args, z = z,
                          lower = -Inf, upper = Inf)$value
  }

  return(aux_fun(x, ..., z = z) - log(z))

}

#### Univariate Tempered Finite Mixture ####

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

