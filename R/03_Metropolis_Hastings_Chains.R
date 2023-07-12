
#' Random Walk Metropolis and MH Chain
#'
#' Functions for running a single chain using the Random Walk Metropolis or a general Metropolis Hastings
#' algorithm for a user-specified target.
#'
#' Missing details for dimension and burn-in
#' # RW Sampler
#'
#' Missing section explaining behavior of `scale` and default/custom samplers.
#'
#' @details
#'  The `_leaner `functions, for memory performance considerations, don't store all the information general ones do.
#'  In particular, they only return the chain values `x` and the global acceptance rate `acc_rate`.
#'
#' @param l_target Function which returns the log-density of the target distribution (possibly up to a constant)
#' and which accepts as first argument the state of the chain.
#' @param ... (Optional) further arguments passed to `l_target`.
#' @param scale Scale of the RWM proposal (see `RW Sampler`)
#' @param S Number of samples of the chain, including burn-in but excluding starting state.
#' @param burn (Default = 0) Number of burn-in iterations. By default, the starting state of the chain
#' is always burned; if you would like to keep it, set this argument to -1 (see `Details`).
#' @param x_0 (Optional) Starting state of the chain.
#' @param x_0_u (Default = 2) If `x_0` is not specified, then the starting state of the chain is randomly sampled
#' from a uniform distribution on (-`x_0_u`,`x_0_u`). If `x_0` is present, it is ignored.
#' @param l_0 (Optional) Precomputed value of the log-density at `x_0`.
#' @param seed (Optional) Seed for reproducibility of the algorithm passed on to `set.seed()`.
#' @param custom_rw_sampler (Optional) User-specified Random Walk sampler function, must accept as its first two arguments
#' the current state of the chain and `scale` (see `RW Sampler`).
#' @param more_sampler_args (Optional) A list of further arguments passed on to `custom_rw_sampler`.
#' @param d (Optional) Dimensionality of the problem (see `Details`).
#' @param silent (Default = FALSE) Should the algorithm avoid printing messages?
#'
#' @return A list containing the results of the Random-Walk Metropolis Chain:
#'
#' * `x`:
#'    A matrix whose rows contain the samples of the chain (after burn-in).
#' * `l_x`:
#'    A vector of log-densities of `x`.
#' * `acc`:
#'    A vector specifying whether or not each of the proposals was accepted or rejected.
#' * `acc_rate`:
#'    Acceptance rate of the chain.
#' * `y`:
#'    A matrix whose rows contain the proposals made during the chain (excluding burn-in).
#' * `l_y`:
#'    A vector of log-densities of `y`.
#' * `delta_l`:
#'    A vector of MH log-ratios.
#'
#' @export
#'
#' @examples
rwm_sampler_chain <- function(l_target, ...,
                              scale = 2.38, S = 1000, burn = 0,
                              x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                              custom_rw_sampler = NULL, more_sampler_args = NULL,
                              d = NULL, silent = FALSE){

#--- Preparation -------------

  # Define the dimension and sampler scale
  stopifnot(is.numeric(scale))
  d <- d %||% ifelse(is.matrix(scale), nrow(scale), length(scale))
  if(d > 1 && !is.matrix(scale)){
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
    scale <- diag(scale)
  }

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # If the user didn't, we define the proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  y <- matrix(nrow = S, ncol = d)
  l_x <- numeric(S + 1)
  l_y <- numeric(S)
  delta_l <- numeric(S)
  acc <- logical(S)

#--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x[1] <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step(x[s, ], l_x[s], l_target, ...,
                                         sampler = sampler,
                                         sampler_args = c(list(scale), more_sampler_args),
                                         do_checks = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_x[s+1] <- rwm_step$l_next
    y[s, ] <- rwm_step$x_prop
    l_y[s] <- rwm_step$l_prop
    delta_l[s] <- rwm_step$delta_l
    acc[s] <- rwm_step$accepted
  }

#--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(mget(c("x", "l_x", "acc", "acc_rate", "y", "l_y", "delta_l")))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x = l_x[-1],
                acc = acc, acc_rate = acc_rate,
                y = y, l_y = l_y, delta_l = delta_l))
  }

  bwx <- 1:(burn+1)
  bwy <- 1:burn
  return(list(x = x[-bwx, ], l_x = l_x[-bwx],
              acc = acc[-bwy], acc_rate = acc_rate,
              y = y[-bwy, ], l_y = l_y[-bwy], delta_l = delta_l[-bwy]))

}

#'
#' @rdname rwm_sampler_chain
#'
#' @param mh_sampler (Optional) User-specified Metropolis-Hastings sampler function,
#' must accept as its first argument the current state.
#' @param other_sampler_args (Optional) A list of further arguments passed on to `mh_sampler`.
#' @param lq_mh (Optional) User-specified log-density function for the MH sampler,
#' must accept as its first argument the current/proposed state.
#' @param other_lq_args (Optional) A list of further arguments passed on to `lq_mh`.
#'
#' @export
#'
mh_sampler_chain <- function(l_target, ...,
                             mh_sampler, other_sampler_args = NULL,
                             lq_mh, other_lq_args = NULL,
                             S = 1000, burn = 0, d = 1,
                             x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, silent = FALSE){

#--- Preparation -------------

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  y <- matrix(nrow = S, ncol = d)
  l_x <- numeric(S + 1)
  l_y <- numeric(S)
  delta_l <- numeric(S)
  acc <- logical(S)

#--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x[1] <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- mh_sampling_step(x_curr = x[s, ], l_curr = l_x[s],
                                 l_target, ...,
                                 sampler = mh_sampler, sampler_args = other_sampler_args,
                                 lq_sampler = lq_mh, lq_sampler_args = other_lq_args,
                                 do_checks = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_x[s+1] <- rwm_step$l_next
    y[s, ] <- rwm_step$x_prop
    l_y[s] <- rwm_step$l_prop
    delta_l[s] <- rwm_step$delta_l
    acc[s] <- rwm_step$accepted
  }

#--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(mget(c("x", "l_x", "acc", "acc_rate", "y", "l_y", "delta_l")))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x = l_x[-1],
                acc = acc, acc_rate = acc_rate,
                y = y, l_y = l_y, delta_l = delta_l))
  }

  bwx <- 1:(burn+1)
  bwy <- 1:burn
  return(list(x = x[-bwx, ], l_x = l_x[-bwx],
              acc = acc[-bwy], acc_rate = acc_rate,
              y = y[-bwy, ], l_y = l_y[-bwy], delta_l = delta_l[-bwy]))

}


#' @rdname rwm_sampler_chain
#'
#' @export
#'
#'
rwm_sampler_leaner_chain <- function(l_target, ...,
                                     scale = 2.38, S = 1000, burn = 0,
                                     x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                                     custom_rw_sampler = NULL, more_sampler_args = NULL,
                                     d = NULL, silent = FALSE){

  #--- Preparation -------------

  # Define the dimension and sampler scale
  stopifnot(is.numeric(scale))
  d <- d %||% ifelse(is.matrix(scale), nrow(scale), length(scale))
  if(d > 1 && !is.matrix(scale)){
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
    scale <- diag(scale)
  }

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # If the user didn't, we define the proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step(x[s, ], l_x, l_target, ...,
                                         sampler = sampler,
                                         sampler_args = c(list(scale), more_sampler_args),
                                         do_checks = FALSE, full_return = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  return(list(x = x[-(1:(burn+1)), ], l_x_curr = l_x, acc = acc[-seq(1,burn)], acc_rate = acc_rate))

}

#' @rdname rwm_sampler_chain
#'
#' @export
#'
#'
mh_sampler_leaner_chain <- function(l_target, ...,
                                    mh_sampler, other_sampler_args = NULL,
                                    lq_mh, other_lq_args = NULL,
                                    S = 1000, burn = 0, d = 1,
                                    x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, silent = FALSE){

  #--- Preparation -------------

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- mh_sampling_step(x_curr = x[s, ], l_curr = l_x,
                                 l_target, ...,
                                 sampler = mh_sampler, sampler_args = other_sampler_args,
                                 lq_sampler = lq_mh, lq_sampler_args = other_lq_args,
                                 do_checks = FALSE, full_return = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  return(list(x = x[-seq(1,burn+1), ], l_x_curr = l_x, acc = acc[-seq(1,burn)], acc_rate = acc_rate))

}

rwm_sampler_leaner_chain_list <- function(l_target, ...,
                                          scale = 2.38, S = 1000, burn = 0,
                                          x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                                          custom_rw_sampler = NULL, more_sampler_args = NULL,
                                          d = NULL, silent = FALSE){

  #--- Preparation -------------

  # Define the dimension and sampler scale
  stopifnot(is.numeric(scale))
  d <- d %||% ifelse(is.matrix(scale), nrow(scale), length(scale))
  if(d > 1 && !is.matrix(scale)){
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
    scale <- diag(scale)
  }

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # If the user didn't, we define the proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step_cpp(x[s, ], l_x, l_target, ...,
                                             sampler = sampler,
                                             sampler_args = c(list(scale), more_sampler_args))
    x[s+1, ] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  return(list(x = x[-(1:(burn+1)), ], l_x_curr = l_x, acc = acc[-seq(1,burn)], acc_rate = acc_rate))

}

mh_sampler_leaner_chain_list <- function(l_target, ...,
                                         mh_sampler, other_sampler_args = NULL,
                                         lq_mh, other_lq_args = NULL,
                                         S = 1000, burn = 0, d = 1,
                                         x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, silent = FALSE){

  #--- Preparation -------------
  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- mh_sampling_step_cpp(x_curr = x[s, ], l_curr = l_x,
                                      l_target, ...,
                                      sampler = mh_sampler, sampler_args = other_sampler_args,
                                      lq_sampler = lq_mh, lq_sampler_args = other_lq_args)
    x[s+1, ] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  return(list(x = x[-seq(1,burn+1), ], l_x_curr = l_x, acc = acc[-seq(1,burn)], acc_rate = acc_rate))

}

rwm_sampler_leaner_chain_list_byprod <- function(l_target_byprod, ...,
                                                 scale = 2.38, S = 1000, burn = 0,
                                                 x_w_0 = NULL, x_0_u = 2,
                                                 l_0 = NULL, by_prod_0 = NA_real_,
                                                 seed = NULL,
                                                 custom_rw_sampler = NULL, more_sampler_args = NULL,
                                                 d = NULL, dim_byprod = 1,
                                                 silent = FALSE){

  #--- Preparation -------------

  # Define the dimension and sampler scale
  stopifnot(is.numeric(scale))
  d <- d %||% ifelse(is.matrix(scale), nrow(scale), length(scale))
  if(d > 1 && !is.matrix(scale)){
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
    scale <- diag(scale)
  }

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # If the user didn't, we define the proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  print(x_w_0)
  if(is.null(x_w_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_w_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_w_0) && length(x_w_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  by_prod <- matrix(nrow = S + 1, ncol = dim_byprod)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_w_0
  if(is.null(l_0)){
    l_x <- l_0
    by_prod[1, ] <- by_prod_0
  }else{
    l_eval <- l_target_byprod(x_0, ...)
    l_x <- l_eval$l_eval
    by_prod[1, ] <- l_eval$by_prod
  }

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step_list_byprod(x[s, ], l_x,
                                                     l_target_byprod, ...,
                                                     sampler = sampler,
                                                     sampler_args = c(list(scale), more_sampler_args))
    x[s+1, ] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
    by_prod[s+1, ] <- rwm_step$by_prod
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate, by_prod = by_prod))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x_curr = l_x, acc = acc, acc_rate = acc_rate, by_prod = by_prod[-1, ]))
  }

  return(list(x = x[-(1:(burn+1)), ], l_x_curr = l_x,
              acc = acc[-seq(1,burn)], acc_rate = acc_rate,
              by_prod = by_prod[-(1:(burn+1)), ]))

}

mh_sampler_leaner_chain_list_byprod <- function(l_target_byprod, ...,
                                                mh_sampler, other_sampler_args = NULL,
                                                lq_mh, other_lq_args = NULL,
                                                S = 1000, burn = 0,
                                                d = 1, dim_byprod = 1,
                                                x_0 = NULL, x_0_u = 2,
                                                l_0 = NULL, by_prod_0 = NA_real_,
                                                seed = NULL, silent = FALSE){

  #--- Preparation -------------
  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  by_prod <- matrix(nrow = S + 1, ncol = dim_byprod)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  if(is.null(l_0)){
    l_x <- l_0
    by_prod[1, ] <- by_prod_0
  }else{
    l_eval <- l_target_byprod(x_0, ...)
    l_x <- l_eval$l_eval
    by_prod[1, ] <- l_eval$by_prod
  }

  # Run iterations
  for(s in 1:S){
    rwm_step <- mh_sampling_step_list_byprod(x_curr = x[s, ], l_curr = l_x,
                                             l_target_byprod, ...,
                                             sampler = mh_sampler, sampler_args = other_sampler_args,
                                             lq_sampler = lq_mh, lq_sampler_args = other_lq_args)
    x[s+1, ] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
    by_prod[s+1, ] <- rwm_step$by_prod
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate, by_prod = by_prod))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x_curr = l_x, acc = acc, acc_rate = acc_rate, by_prod = by_prod[-1, ]))
  }

  return(list(x = x[-seq(1,burn+1), ], l_x_curr = l_x,
              acc = acc[-seq(1,burn)], acc_rate = acc_rate,
              by_prod = by_prod[-seq(1,burn+1), ]))

}

