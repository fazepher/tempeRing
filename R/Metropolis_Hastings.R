
checks_mh <- function(x_curr, x_prop, l_curr, l_prop, lq_c2p, lq_p2c, C, d, x_is_matrix){

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


  return(TRUE)

}

#' Metropolis-Hastings Step
#'
#' @param x_curr A vector or matrix of current state(s)
#' @param x_prop A vector or matrix of proposed state(s)
#' @param l_curr A vector of current log-densities, possibly up to proportionality
#' @param l_prop A vector of proposed log-densities, possibly up to proportionality
#' @param lq_c2p A vector of transition log-densities from the current state(s) to the proposal(s)
#' @param lq_p2c A vector of transition log-densities from the porposal(s) to the currrent state(s)
#' @param do_checks If TRUE (the default), run preliminary checks for validity of arguments
#'
#' @return A list containing the result of the Metropolis-Hastings step
#' @export
#'
mh_step <- function(x_curr, x_prop, l_curr, l_prop, lq_c2p = 0, lq_p2c = 0, do_checks = TRUE){

  # Preparation
  ret_list <- c("x_next", "l_next", "accepted")
  C <- length(l_curr)
  x_is_matrix <- is.matrix(x_curr)
  d <- ifelse(x_is_matrix, ncol(x_curr),
              ifelse(C==1, length(x_curr), 1))

  # Checks if requested
  if(do_checks){
    checks_good <- checks_mh(x_curr, x_prop, l_curr, l_prop, lq_c2p, lq_p2c, C, d, x_is_matrix)
  }

  # Acceptance
  delta_l <- (l_prop + lq_p2c) - (l_curr + lq_c2p)
  accepted <- vapply(X = delta_l,
                     FUN = function(w){ w > 0 || log(stats::runif(1)) <= w },
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

metropolis_step <- function(x_curr, x_prop, l_curr, l_prop, ...){

  mh_step(x_curr,x_prop,l_curr,l_prop, ...)

}
