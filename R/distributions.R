
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
#' ## Vectorization
#'
#' The functions inherit vectorization, so can be used for example with
#' `ggplot2::geom_function()` to show the flattening effect of tempering
#'
#' ```
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
#' ```
#'
#'
#' @param x Vector of quantiles.
#' @param mean Vector of means.
#' @param sd Vector of standard deviations.
#' @param beta Inverse temperature parameter \eqn{\beta > 0}.
#' @param n Number of observations.
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
#' @inheritParams rnorm_temp
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

#' General Finite Mixture log densities and samplers
#'
#' Log-density, density, and random generation of a mixture of a user specified
#' common parametric **log-density** (or distribution for sampling), according to some
#' mixture weights `w`.
#'
#' In general, suppose we have a parametric family of distributions such that \eqn{f(x | \theta)}
#' is the density given the parameter vector \eqn{\theta}. A K component mixture of this
#' distribution with non-negative weights \eqn{w_1,...,w_K} such thath \eqn{\sum w_k = 1} has
#' density
#' \deqn{f(x|w,\Theta) = \sum w_k f(x|\theta_k)}
#' where \eqn{\theta_k} represent the specific parameters that the k-th component takes.
#'
#' For numerical reasons, it is recommended to work on the log-scale so that for a desired
#' distribution, the corresponding mixture can be generated from the log-densities
#' \eqn{l(x|\theta) = log( f(x | theta))}, via de Log-Sum-Exp form (see [matrixStats::logSumExp]):
#' \deqn{l(x|w,\Theta) = LSE[log(w_k) + l(x)|\theta_k)] over k = 1,... K}
#'
#' # Component Parameters
#'
#' Each of the mixture components would have their specific values of the parameters \eqn{\theta_k}.
#' But some parameters may remain fixed across components. For example, suppose we want to work
#' with a mixture of normal distributions in which all components have the same standard deviation
#' but different means. There are, at least, two ways of evaluating a log-density in this scenario;
#' namely, repeat the shared parameters and treat them as if they were varying through `...`
#' or provide them via the `shared_args` list:
#'
#' ```
#' lmix(x = 0, w = c(0.5, 0.5), ldens = lnorm,
#'      mean = c(-1, 1), sd = c(1, 1))
#'
#' lmix(x = 0, w = c(0.5, 0.5), ldens = lnorm,
#'      mean = c(-1,1), shared_args = list(sd = 1))
#'
#' ```
#'
#' Component parameters fed through `...` are passed on to `mapply`, so they need to be
#' unambigously iterable across mixture components, thus it is recommended to pass them as lists
#' if we want to avoid surprises. For instance, if now we wanted to evaluate a density for
#' a mixture of multivariate normals,  we can pass the means and covariances as lists of k
#' elements each.
#'
#' ```
#' dmix(x=c(0,0), w = c(0.2, 0.8), ldens = lmvtnorm,
#'      mu = list(c(-3, -3), c(1, 1)),
#'      sigma = list(diag(0.5,2), diag(0.5,2)))
#'
#' dmix(x=c(0,0), w = c(0.2, 0.8), ldens = lmvtnorm,
#'      mu = list(c(-3,3), c(1,1)),
#'      shared_args = list(sigma = diag(0.5,2)))
#'
#' ```
#'
#' The behavior for generating random samples is different, though. As each sample is only
#' generated from a given component and there is, in general, no need to sweep through the parameters
#' of each component. For this reason, in this case component parameters are passed as a single linst
#' in which each element contains the parameters of each component. In the multivariate normal mixture
#' scenario this is:
#'
#' ```
#' rmix(n=3, w = c(0.2, 0.8), rdist = rmvtnorm,
#'      comp_args_list = list(k1 = list(mu = c(-3, 3)),
#'                            k2 = list(mu = c(1, 1))),
#'      shared_args = list(sigma = diag(0.5,2)))
#'
#' rmix(n=3, w = c(0.2, 0.8), rdist = rmvtnorm,
#'      comp_args_list = list(k1 = list(mu = c(-3, 3)),
#'                            k2 = list(mu = c(1, 1))),
#'      shared_args = list(sigma = diag(0.5,2)),
#'      simplify = TRUE)
#'
#' ```
#'
#'
#' @param x A quantile or quantiles if `ldens` allows for vectorization.
#' @param w A vector of non-negative mixture weights. To be a valid mixture they must sum to 1.
#' @param ldens A function that returns the log-density of the desired common mixture distribution.
#' @param ... Other arguments passed on to `ldens` that vary by mixture component
#' (see Component Parameters).
#' @param shared_args List of other arguments passed on to `ldens` who are shared
#' by all mixture components (see Component Parameters).
#'
#' @seealso [lmix_norm], [lmix_mvtnorm], [lmix_skewnorm] for some specific mixtures.
#' For some tempered mixtures see [ulmix_norm_temp] or [ulmix_temp] in general.
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
#' @param log For `dmix`, whether to return the log-density or not (the default).
#'
#' @export
#'
dmix <- function(x, w, ldens, ..., shared_args = NULL, log = FALSE){
  l <- lmix(x, w, ldens, ..., shared_args = shared_args)
  if(log) return(l) else return(exp(l))
}
#' @rdname lmix
#'
#' @inheritParams rnorm_temp
#' @param rdist For `rmix`, sampler function of the desired parametric distribution, whose first
#' argument must be named `n` and refer to a number of samples asked from it.
#' @param comp_args_list For `rmix`, a list whose k-th element contains the
#' list of parameters specific to the k-th mixture component (see Component Parameters).
#' @param simplify For `rmix` whether to attempt to simplify the list of
#' generated samples into an array or not (the default).
#'
#' @export
#'
rmix <- function(n, w, rdist, comp_args_list, shared_args = NULL, simplify = FALSE){

  # Select components according to w
  k_vec <- sample.int(length(w), size = n, replace = TRUE, prob = w)

  # Sample from the distribution
  if(simplify){
    x <- sapply(k_vec, function(k) do.call(rdist, c(list(n=1), comp_args_list[[k]], shared_args)),
                simplify = TRUE)
  }else{
    x <- lapply(k_vec, function(k) do.call(rdist, c(list(n=1), comp_args_list[[k]], shared_args)))
  }

  return(x)
}

#### Finite Mixture of Normals ####

#' Finite Mixture of Normal distributions
#'
#' Log-density, density, random sampling from a finite mixture of normal distributions,
#' as well as unnormalized log-density of the tempered mixture.
#'
#' The tempering in `ulmix_norm_temp` is done *after* the mixture; this is not to be confused
#' with a mixture of tempered normals. This latter scenario would just be a regular mixture of normals,
#' corresponding to `lmix_norm`, as tempering a normal just rescales it (see [lnorm_temp]).
#'
#'
#' @inheritParams lnorm
#' @inheritParams lmix
#' @param mean Vector containing the means of the components
#' @param sd (Optional) Vector containing the standard deviations of the components
#' @param shared_sd Shared standard deviation of the components; for different scales use `sd`.
#'
#' @returns
#'
#' `lmix_norm` gives the log-density, `dmix_norm` the density and `rmix_norm` generates random samples.
#' Finally, `ulmix_norm_temp` gives the unnormalized log-density of the tempered mixture.
#'
#' @export
#'
lmix_norm <- function(x, w, mean, sd = NULL, shared_sd = 1){
  if(is.null(sd)){
    lmix(x, w, lnorm, mean, shared_args = list(sd = shared_sd))
  }else{
    lmix(x, w, lnorm, mean, sd)
  }
}

#' @rdname lmix_norm
#'
#' @param log For `dmix_norm`, whether to return the log-density or not (the default).
#'
#' @export
#'
dmix_norm <- function(x, w, mean, sd = NULL, shared_sd = 1, log = FALSE){
  if(is.null(sd)){
    dmix(x, w, lnorm, mean, shared_args = list(sd = shared_sd), log = log)
  }else{
    dmix(x, w, lnorm, mean, sd, log = log)
  }
}

#' @rdname lmix_norm
#'
#' @inheritParams rnorm_temp
#'
#' @export
#'
rmix_norm <- function(n, w, mean, sd = NULL, shared_sd = 1){
  if(is.null(sd)){
    rmix(n, w, rnorm, mean, shared_args = list(sd = shared_sd))
  }else{
    rmix(n, w, rnorm, mean, sd)
  }
}

#' @rdname lmix_norm
#'
#' @inheritParams lnorm_temp
#'
#' @export
#'
ulmix_norm_temp <- function(x, beta = 1, w, mean, sd = NULL, shared_sd = 1){
  if(is.null(sd)){
    beta*lmix(x, w, lnorm, mean, shared_args = list(sd = shared_sd))
  }else{
    beta*lmix(x, w, lnorm, mean, sd)
  }
}

#### Finite Mixture of Multivariate Normals ####

lmix_mvtnorm <- function(x, w, ..., shared_args = NULL){
  lmix(x, w, lmvtnorm, ..., shared_args = shared_args)
}
dmix_mvtnorm <- function(x, w, ..., shared_args = NULL, log = FALSE){
  dmix(x, w, dmvtnorm, ..., shared_args = shared_args, log = log)
}
rmix_mvtnorm <- function(n, w, comp_args_list, shared_args = NULL, simplify = FALSE){
  rmix(n, w, rdist = rmvtnorm,
       comp_args_list = comp_args_list, shared_args = shared_args,
       simplify = simplify)
}
ulmix_mvtnorm_temp <- function(x, beta = 1, w, ..., shared_args = NULL){
  beta*lmix_mvtnorm(x, w, ..., shared_args = shared_args)
}

#### General tempered densities ####

# Mostly for completeness and as a simpler version to code than the mixture ones,
# hence they are not exported for now.
ul_temp <- function(x, beta = 1, ldens, ...){
  beta*ldens(x, ...)
}
l_temp <- function(x, beta = 1, ldens, ..., z = NULL){

  aux_fun <- function(y, aux_beta = beta, aux_ldens = ldens, ...,
                      aux_dens_is_log = dens_is_log, lz = lz){

    ul <- ul_temp(y, aux_beta, aux_ldens, ..., dens_is_log = aux_dens_is_log)

    if(is.null(z)) return(exp(ul)) else return(ul)

  }

  if(is.null(z)){
    z <- stats::integrate(aux_fun, aux_beta = beta, aux_ldens = ldens, ..., z = z,
                          lower = -Inf, upper = Inf)$value
  }

  return(aux_fun(x, ..., z = z) - log(z))

}

#### General Tempered Mixtures ####

#'
#' Tempered Mixture log densities generator
#'
#' @description
#' Unnormalized and normalized log-tempered densities from a mixture of a user defined
#' **log-density**.
#'
#' While mixing densities still yields a normalized density, tempering *does not* do so
#' and requires renormalization. Hence, the direct result of tempering is the unnormalized
#' `ulmix_temp`. The normalized version `lmix_temp` is provided mainly for **univariate** log-densities,
#' where internal numerical integration is used via `stats::integrate`.
#'
#' ## Component Parameters
#'
#' Same behavior as [lmix]: `...` expects iterable arguments to be swept across all components while
#' `shared_args` is a list of shared arguments across components.
#' For a more detailed explanation please see the documentation and examples in [lmix].
#'
#' @details
#' When using `lmix_temp` for general *multivariate* mixtures, the user should provide the
#' normalizing constant `log_z` (in the log scale). This parameter can also be useful for efficiency in
#' one-dimensional situations, as one would avoid estimating it at every call.
#' Note, however, that many tasks in the context of the package don't necessarily require the
#' normalized version; this function is provided for completeness and as a tool
#' for those situations were the normalized versions are indeed sought.
#'
#' @inheritParams lmix
#' @param beta Inverse temperature parameter β > 0.
#' @param ldens A function that returns the log-density of the desired common mixture distribution.
#' For `lmix_temp` see Details for restrictions on multivariate densities.
#'
#'
#' @returns
#' `ulmix_temp` returns the unnormalized log-density of the mixture, while `lmix_temp` estimates
#' the normalized value.
#'
#' @seealso [lmix] [ulmix_norm_temp]
#'
#' @export
#'
ulmix_temp <- function(x, beta = 1, w, ldens, ..., shared_args = NULL){
  beta*lmix(x, w, ldens, ..., shared_args = shared_args)
}

#' @rdname ulmix_temp
#'
#' @param log_z For `lmix_temp`, normalizing constant. If NULL (the default) it is estimated
#' via `stats::integrate`. Hence, while optional for univariate mixtures, it is necessary
#' for proper behavior on multivariate mixtures (see Details).
#'
#' @export
#'
lmix_temp <- function(x, beta = 1, w, ldens, ..., shared_args = NULL, log_z = NULL){

  aux_fun <- function(x, aux_beta = beta, aux_w = w,  aux_ldens = ldens, ...,
                      aux_shared_args = shared_args, log_z = log_z){

    ul <- ulmix_temp(x, aux_beta, aux_w, aux_ldens, ..., shared_args = aux_shared_args)

    if(is.null(log_z)) return(exp(ul)) else return(ul)

  }

  if(is.null(log_z)){
    log_z <- stats::integrate(aux_fun, aux_beta = beta, aux_w = w,
                              aux_ldens = ldens, ..., aux_shared_args = shared_args,
                              log_z = log_z, lower = -Inf, upper = Inf)$value |>
      log()
  }

  return(aux_fun(x, ..., log_z = log_z) - log_z)

}

