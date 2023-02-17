
#' @export
get_HAT_info <- function(mode_guess, l_target, ..., beta_hat = NULL,
                         optimize = TRUE, method = "BFGS", control_optim = list(fnscale = -1)){

  bh <- beta_hat %||% 1
  if(optimize){
    optimizations <- lapply(mode_guess,
                            function(m) optim(m, l_target, beta = bh, ..., method = method,
                                              control = control_optim, hessian = TRUE))
    lapply(optimizations, function(o) c(o$convergence, o$counts)) |> print()
    modes <- lapply(optimizations,function(o) o$par)
    mH <- lapply(optimizations, function(o) -o$hessian)
  }else{
    modes <- mode_guess
    mH <- lapply(modes, function(m) -(optimHess(m, l_target, ...)))
  }

  l_target_modes <- vapply(modes, l_target, numeric(1), beta = 1, ...)

  cholCov_inv <- lapply(mH, chol)
  hldetCov_inv <- vapply(cholCov_inv, function(chS) sum(log(diag(chS))), numeric(1))
  Cov <- lapply(mH,solve)
  cholCov <- lapply(Cov,chol)
  half_l_detCov <- vapply(cholCov, function(chS) sum(log(diag(chS))), numeric(1))
  detCov <- exp(2*half_l_detCov)

  l_w_tilde <- l_target_modes + half_l_detCov
  ls_w_tilde <- matrixStats::logSumExp(l_w_tilde)
  w <- exp(l_w_tilde-ls_w_tilde)

  return(mget(c("mH","cholCov_inv","hldetCov_inv",
                "Cov","cholCov","half_l_detCov","detCov",
                "l_target_modes","w","modes")))

}

modAssignment <- function(x, beta, HAT_info, assign_type){
  if(assign_type == "euclidean"){
    modAssignment_euclidean(x, beta, HAT_info)
  }
}

modAssignment_euclidean <- function(x, beta, HAT_info){

  lP_vec <-mapply(function(mu,w,Cov_inv, hlds) log(w) + beta*lmvtnorm(x, mu, sigma_inv = Cov_inv, logdet_sigma = 2*hlds),
                  mu = HAT_info$modes, w = HAT_info$w, Cov_inv = HAT_info$mH, hlds = HAT_info$hldetCov_inv)
  A <- which.max(lP_vec)

  return(list("A" = A[1], "lP_j" = lP_vec[A[1]]))

}

#' @export
lHAT_target <- function(x, beta, HAT_info, ltemp_target, ..., silent = FALSE){

  ## Basic Weight Preservation

  # Early exit if beta is 1
  l_eval <- ltemp_target(x, beta = beta, ...)
  if(beta == 1){ return(l_eval) }

  # Assign the mode at beta and 1
  mod_beta <- modAssignment_euclidean(x, beta, HAT_info)
  mod_1 <- modAssignment_euclidean(x, 1, HAT_info)
  l_mod <- HAT_info$l_target_modes[ mod_beta$A ]
  if(is.na(mod_beta$lP_j) || is.na(mod_1$lP_j) || is.na(l_mod)){stop("Error")}

  # Early exit if the assigned modes are the same
  if(mod_beta$A == mod_1$A){
    return(l_eval + (1-beta)*l_mod)
  }

  ## G Weight Preservation only when the assigned modes differ
  l_corr_ctes <- 0.5*length(HAT_info$modes[[1]])*(log(2*pi) - log(beta))
  l_approx_w <- mod_beta$lP_j - log(HAT_info$w[mod_beta$A])
  l_G <- l_mod + l_corr_ctes + 1/HAT_info$hldetCov_inv[mod_beta$A] + l_approx_w
  return(l_G)

}

#' @export
HAT_rwm_chain <- function(ltemp_target, ..., HAT_info,
                          beta_schedule, swap_type = "deo",
                          scale = 1, Cycles = 1000, Temp_Moves = 1, Within_Moves = 5,
                          x_0 = NULL, x_0_u = 2, seed = NULL,
                          custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                          silent = FALSE){

  # Then it's just a wrapper of the PT algorithm for the new HAT target
  hat_args <- c(list(l_target = lHAT_target, HAT_info = HAT_info,
                     ltemp_target = ltemp_target), rlang::dots_list(...),
                list(beta_schedule = beta_schedule, swap_type = swap_type,
                     scale = scale, Cycles = Cycles, Temp_Moves = Temp_Moves, Within_Moves = Within_Moves,
                     x_0 = x_0, x_0_u = x_0_u, seed = seed,
                     custom_rw_sampler = custom_rw_sampler, target_names = target_names, d = d,
                     silent = silent))
  do.call(PT_rwm_chain, hat_args)

}
