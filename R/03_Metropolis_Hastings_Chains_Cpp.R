
#' Random-Walk Metropolis
#'
#' A single global scale real-valued parameter is assumed, dimension must be explicitly passed.
#' Less checks as this is mainly for developing ALPS internally
#'
#' @param l_target Function which returns the log-density of the target distribution (possibly up to a constant)
#' and which accepts as first argument the state of the chain.
#' @param d Integer specifying the dimensionality of the state vector.
#' @param ... Additional parameters passed to `l_target`.
#' @param global_scale Scale for the RWM sampler to use.
#' @param S Number of total samples to run the Markov Chain for.
#' @param burn Number of simulated samples to discard as burn-in. If `0` (the default),
#' the initial state `x_0` is discarded but otherwise all the simulated samples are kept.
#' Use `-1` to also keep he initial state `x_0`.
#' @param x_0 (Optional) Initial state of the chain.
#' @param x_0_u (Optional) Positive number indicating the upper limit of a Uniform distribution
#' used to sample an initial value for the chain. The lower limit is then its negative. Hence,
#' each component of `x_0` is sampled uniformly between -`x_0_u` and `x_0_u`.
#' Ignored if `x_0` is present.
#' @param l_0 (Optional) Precomputed value of `l_target` evaluated at `x_0`.
#' @param seed (Optional) Seed to set before starting sampling the chain.
#' @param custom_rw_sampler (Optional) Function for a user defined Random-Walk sampler.
#' It must accept as its first argument the current state of the chain and as its second the `global_scale`
#' used for the chain. If `NULL` (the default) the chain assumes a Normal RW sampler, see `Details`.
#' @param more_sampler_args (Optional) List containing further arguments passed to `custom_rw_sampler`.
#' @param silent (Logical) Should message-printing be disabled? By default, it is `FALSE`.
#'
#' @returns A list containing the results of the RWM Chain:
#' list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate)
#' * `x`:
#'    A matrix containing the evolution of the states of the chain with `d` rows and `S` - `burn`
#'    columns.
#' * `l_x_curr`:
#'    Current (last) value of the log-density of the target.
#' * `acc`:
#'    A vector specifying whether or not each proposal was accepted or rejected (excluding burn-in).
#' * `acc_rate`:
#'    Observed acceptance rate across all iterations (including burn-in).
#'
rwm_global_scale_sampler_leaner_chain <- function(l_target, d, ...,
                                                  global_scale = 2.38, S = 1000, burn = 0,
                                                  x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                                                  custom_rw_sampler = NULL, more_sampler_args = NULL,
                                                  silent = FALSE){

  #--- Preparation -------------

  # If the user didn't, we define the proposal sampler as indep. normals
  if(is.null(custom_rw_sampler)){
    sampler <- function(x,scale){
      rmvtnorm_temp_cpp(n = 1, mu = x, global_scale = scale)
    }
  }else{
    sampler <- custom_rw_sampler
  }

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }

  # Preallocate containers
  x <- matrix(nrow = d, ncol = S + 1)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[, 1] <- x_0
  l_x <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step_cpp(x[, s], l_x, l_target, ...,
                                             sampler = sampler,
                                             sampler_args = c(list(global_scale), more_sampler_args))
    x[, s+1] <- rwm_step$x_next
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
    return(list(x = x[,-1], l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  return(list(x = x[,-(1:(burn+1))], l_x_curr = l_x, acc = acc[-seq(1,burn)], acc_rate = acc_rate))

}
