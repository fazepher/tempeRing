
#' Metropolis-Hastings Step
#'
#' Perform a Metropolis-Hastings step on one or more states.
#' * If we already have proposed states and all log-density factors we can call `mh_step()`.
#' * If we wish to first sample a proposal and then give the step, we can use `mh_sampling_step()`.
#' * Both `metropolis_step()` and `metropolis_sampling_step()` are just wrappers for symmetric
#' transition kernels.
#'
#' # Dimension
#'
#' We want to perform MH steps on `C` states, so all log-density inputs are vectors of length `C`.
#' The only exception is for the transition log-densities that by default are set both to `0`
#' to signal that we are dealing with symmetrical Metropolis steps.
#' Similarly, the `accepted` and `l_next` elements of the output list are always vectors of length `C`.
#'
#' States, on the other hand, live in a space of dimension `d`. There are thus 3 big combinations:
#'
#' * `C`=1: We want to perform a step on a single d-dimensional state. Without loss of generality,
#' `x_curr` and `x_prop` are two vectors of length `d`.
#' * `C`>1, `d`>1: We want to perform a step on each of several multidimensional states.
#' Here, `x_curr` and `x_prop` are two `C` by `d` matrices, so that each state is a row vector.
#' * `C`>1, `d`=1: We want to perform a step on each of several 1-dimensional states.
#' Here, all inputs can be vectors of length `C` or, alternatively,
#' `x_curr` and `x_prop` may be two `C` by `1` matrices.
#'
#' The output element `x_next` is always of the same structure as the input `x_curr`.
#'
#' # The Metropolis-Hastings Ratio
#' The usual form of the MH acceptance probability \eqn{\alpha = min{1, MH ratio}},
#' relies on the ratio
#' \deqn{ MH ratio = \pi(x_1) q(x_0|x_1) / \pi(x_0) q(x_1|x_0) }
#' to satisfy detail balance.
#' For numerical reasons, we wish to work on the log scale and the ratio becomes
#' \deqn{ MH log-ratio =  l(x_1) + lq(x_0|x_1) - l(x_0) - lq(x_1|x_0) }
#' Whenever the transition kernel is symmetrical (i.e. \eqn{q(x_0|x_1)=q(x_1|x_0)})
#' we can omit those terms from the calculation and recover the original Metropolis et. al ratio.
#' This is the default assumption of the `mh_step()` function.
#' ## Note on proportionality
#' We want to also take advantage of the cancellation of normalizing constants so that we only need
#' the log-densities *up to a constant of proportionality*. Note however that this refers
#' to the **same** constant. While there are usually no problems when it comes to the target density,
#' care must be taken with general transition kernels as they may have different
#' underlying densities whose normalizing constants don't cancel.
#'
#'
#' @param x_curr,x_prop Two vectors or matrices of current and proposed states. See `Dimension` for
#'   how to specify these to correctly identify the dimension of the state space.
#' @param l_curr,l_prop Vectors of values of log-densities, possibly up to proportionality, of the
#'   current and proposed states.
#' @param lq_c2p,lq_p2c (Optional) Vectors of transition log-densities, possibly up to proportionality.
#'   The suffix `_c2p` corresponds to transitions from the **c**urrent towards the **p**roposed state(s),
#'   while the suffix `_p2c` contains the reversed transitions. If these values are ommited (the default),
#'   the Hastings ratio becomes standard Metropolis.
#' @param do_checks If `TRUE` (the default), runs some preliminary checks for arguments validity.
#'    It may be set to `FALSE` if one is sure that the arguments are correct and doesn't want to incurr in
#'    slight overhead.
#' @param full_return If `TRUE` (the default), the call returns its full set of results. See  `Value`.
#'
#' @returns A list containing the results of the Metropolis-Hastings step:
#' * `x_next`:
#'    The next values of the chain. This object has the same structure as the input `x_curr`.
#'    For each proposal, it contains the corresponding `x_prop` or `x_curr` values depending on
#'    their acceptance or rejection.
#' * `l_next`:
#'    A vector of log-density values corresponding to `x_next` elements.
#'    These may be re-used in the next iteration of a chain to avoid wasting resources in what
#'    usually is an expensive computation.
#' * `accepted`:
#'    A vector specifying whether or not each of the proposals was accepted or rejected.
#'    May be useful for acceptance rate monitoring.
#' If `full_return` is `TRUE`, then the resulting list also contains the following:
#' * `delta_l`:
#'    The MH log-ratio
#' * The called parameters `x_curr`, `x_prop`, `l_curr`,`l_prop`,`lq_c2p`,`lq_p2c`.
#'
#' @export
#'
mh_step <- function(x_curr, x_prop, l_curr, l_prop, lq_c2p = 0, lq_p2c = 0,
                    do_checks = TRUE, full_return = TRUE){

#--- Preparation and Checks -------------

  # Preparation
  if(full_return){
    ret_list <- c("x_next", "l_next", "accepted", "delta_l",
                  "x_curr", "x_prop", "l_curr", "l_prop", "lq_c2p", "lq_p2c")
  }else{
    ret_list <- c("x_next", "l_next", "accepted")
  }

  C <- length(l_curr)
  x_is_matrix <- is.matrix(x_curr)
  d <- ifelse(x_is_matrix, ncol(x_curr),
              ifelse(C==1, length(x_curr), 1))

  # Checks if requested
  if(do_checks){

    # All numeric
    stopifnot(is.numeric(l_curr))
    stopifnot(is.numeric(l_prop))
    stopifnot(is.numeric(lq_c2p))
    stopifnot(is.numeric(lq_p2c))

    # Valid Lengths
    if(C != length(l_prop)){
      stop("Lengths of l_curr and l_prop must be equal")
    }

    check_l <- length(lq_c2p)
    if(check_l != length(lq_p2c)){
      stop("Lengths of lq_c2p and lq_p2c must be equal")
    }
    if(!check_l %in% c(1,C)){
      stop("lengths of lq_c2p and lq_p2c must match that of l_curr or be 1")
    }

    # X agreements
    if(x_is_matrix){

      if(d != ncol(x_prop)){
        stop("At x_prop, number of columns must match those of x_curr")
      }
      if(C != nrow(x_curr)){
        stop("At x_curr, number of rows must match length of l_curr")
      }
      if(C != nrow(x_prop)){
        stop("At x_prop, number of rows must match length of l_curr")
      }

    } else{
      stopifnot(length(x_curr) == length(x_prop))
    }

  }

#--- Algorithm, Post processing and results -------------

  # Acceptance
  delta_l <- (l_prop + lq_p2c) - (l_curr + lq_c2p)
  accepted <- vapply(X = delta_l,
                     FUN = function(w){ w > 0 || log(runif(1)) <= w },
                     FUN.VALUE = logical(1))

  l_next <- numeric(C)
  l_next[accepted] <- l_prop[accepted]
  l_next[!accepted] <- l_curr[!accepted]


  ## 3 possible return structures for x_next depending on that of x_curr ##

  # X matrix of C d-dimensional states
  if(x_is_matrix){
    x_next <- matrix(nrow = C, ncol = d)
    x_next[accepted, ] <- x_prop[accepted, ]
    x_next[!accepted, ] <- x_curr[!accepted, ]
    return(mget(ret_list))
  }

  # X d-dimensional state single vector
  if(d > 1){
    x_next <- if(accepted){x_prop}else{x_curr}
    return(mget(ret_list))
  }

  # X vector of C 1-dimensional states
  x_next <- numeric(C)
  x_next[accepted] <- x_prop[accepted]
  x_next[!accepted] <- x_curr[!accepted]
  return(mget(ret_list))

}

#' @rdname mh_step
#'
#' @export
#'
metropolis_step <- function(x_curr, x_prop, l_curr, l_prop, do_checks = TRUE, full_return = TRUE){

  mh_step(x_curr, x_prop, l_curr, l_prop, do_checks = do_checks)

}

#' Sample and perform a Metropolis-Hasting Step
#'
#' @rdname mh_step
#'
#' @param l_target For `mh_sampling_step()`, a function that computes the log-density of the target,
#'    possibly up to proportionality. It must accept `x_curr` as its first argument and its returned
#'    value must be a valid corresponding `l_prop`.
#' @param ... Arguments to be passed on to the `l_target` function.
#' @param sampler For `mh_sampling_step()`, a sampling function whose first argument is compatible with
#'    `x_curr` and returns appropriate `x_prop` states.
#' @param sampler_args A list of further arguments to be passed on to the `sampler`
#' @param lq_sampler (Optional) For `mh_sampling_step()`, a function that computes the transition log-densities
#'    for the `sampler()`. It first two arguments must be conformable with `x_curr` and `x_prop` as
#'    it will be used to compute the relevant `lq_c2p` and `lq_p2c` factors of the MH step.
#'    If ommited, it is assumed that it is symmetric and, hence, it's not needed.
#' @param lq_sampler_args A list of further arguments to be passed on to `lq_sampler`.
#'
#' @export
#'
mh_sampling_step <- function(x_curr, l_curr, l_target, ..., sampler, sampler_args = NULL,
                             lq_sampler = NULL, lq_sampler_args = NULL, do_checks = TRUE, full_return = TRUE){

  x_prop <- do.call(sampler, c(list(x_curr), sampler_args))
  l_prop <- do.call(l_target, c(list(x_prop), rlang::dots_list(...)))

  if(is.function(lq_sampler)){
    lq_c2p <- do.call(lq_sampler, c(list(x_curr, x_prop), lq_sampler_args))
    lq_p2c <- do.call(lq_sampler, c(list(x_prop, x_curr), lq_sampler_args))
  }else{
    lq_c2p <- 0
    lq_p2c <- 0
  }

  results <- mh_step(x_curr, x_prop, l_curr, l_prop, lq_c2p, lq_p2c, do_checks = do_checks)

  return(results)

}

#' @rdname mh_step
#'
#' @export
#'
metropolis_sampling_step <- function(x_curr, l_curr, l_target, ..., sampler, sampler_args,
                                     do_checks = TRUE, full_return = TRUE){

  mh_sampling_step(x_curr, l_curr, l_target, ...,
                   sampler = sampler, sampler_args = sampler_args,
                   do_checks = do_checks)

}
