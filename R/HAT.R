
#' @export
get_HAT_info <- function(mode_guess, l_target, ...,
                         optimize = TRUE, method = "Nelder-Mead", control_optim = list(fnscale = -1)){

  if(optimize){
    optimizations <- lapply(mode_guess,
                            function(m) optim(m, l_target, beta = 1, ..., method = method,
                                              control = control_optim, hessian = TRUE))
    modes <- lapply(optimizations,function(o) o$par)
  }else{
    modes <- mode_guess
  }
  mH <- lapply(modes, function(m) -(numDeriv::hessian(l_target, m, ...)))
  Cov <- lapply(mH,solve)
  cholCov <- lapply(Cov,chol)
  detCov <- vapply(cholCov, function(chS) prod(diag(chS))^2, numeric(1))
  l_target_modes <- vapply(modes, l_target, numeric(1), beta = 1, ...)
  w_tilde <- exp(l_target_modes + 0.5*log(detCov))
  w <- w_tilde/sum(w_tilde)

  return(mget(c("mH","Cov","cholCov","detCov","l_target_modes","w","modes")))

}


modAssignment <- function(x, beta, HAT_info){

  lP_vec <-mapply(function(mu,w,Cov) log(w) + lmvtnorm(x, mu, Cov/beta),
                  mu = HAT_info$modes, w = HAT_info$w, Cov = HAT_info$Cov)
  A <- which.max(lP_vec)

  return(list("A" = A[1], "lP_j" = lP_vec[A[1]]))

}

#' @export
lHAT_target <- function(x, beta, HAT_info, ltemp_target, ..., G_type = 1, silent = FALSE){

  ## Basic Weight Preservation

  # Early exit if beta is 1
  l_eval <- ltemp_target(x, beta = beta, ...)
  if(beta == 1){ return(l_eval) }

  # Assign the mode at beta and 1
  mod_beta <- modAssignment(x, beta, HAT_info)
  mod_1 <- modAssignment(x, 1, HAT_info)
  l_mod <- HAT_info$l_target_modes[ mod_beta$A ]
  if(is.na(mod_beta) || is.na(mod_1) || is.na(l_mod)){stop("Error")}

  # Early exit if the assigned modes are the same
  l_w <- l_eval + (1-beta) * l_mod
  if(mod_beta$A == mod_1$A){ return(l_w) }

  ## G Weight Preservation only when the assigned modes differ

  # This appears to be the "non-robust" G, earlier exit if G_type is not robust (the default 1)
  if(G_type != 1){
    return(0.5*length(x)*(log(2*pi)-log(beta)) + mod_beta$lP_j)
  }

  # Robust G otherwise

  weight_G <- (2 / (1 + exp(mod_1$lP_j - mod_beta$lP_j))) - 1

  # Early exit if weight_G is negative
  # I DON'T KNOW THE THEORY HERE, PAPER DOESN'T MENTION IT, ASK NICK AND GARETH
  # I WILL RETURN A -Inf log density to always reject a move to such a state
  if(weight_G < 0){
    if(!silent){warning("weight_G is negative")}
    return(-Inf)
  }

  l_cold <- log(1 - weight_G) + do.call(ltemp_target, list(x=x, beta =1, ...))
  l_wp_beta <- log(weight_G) + l_w
  l_G <- logSumExp(c(l_cold, l_wp_beta))

  return(l_G)

}

#' @export
HAT_rwm_chain <- function(ltemp_target, ..., HAT_info, G_type = 1,
                          beta_schedule, swap_type = "deo",
                          scale = 1, Cycles = 1000, Temp_Moves = 1, Within_Moves = 5,
                          x_0 = NULL, x_0_u = 2, seed = NULL,
                          custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                          silent = FALSE){

  # Then it's just a wrapper of the PT algorithm for the new HAT target
  hat_args <- c(list(l_target = lHAT_target, G_type = G_type, HAT_info = HAT_info,
                     ltemp_target = ltemp_target), rlang::dots_list(...),
                list(beta_schedule = beta_schedule, swap_type = swap_type,
                     scale = scale, Cycles = Cycles, Temp_Moves = Temp_Moves, Within_Moves = Within_Moves,
                     x_0 = x_0, x_0_u = x_0_u, seed = seed,
                     custom_rw_sampler = custom_rw_sampler, target_names = target_names, d = d,
                     silent = silent))
  do.call(PT_rwm_chain, hat_args)

}
