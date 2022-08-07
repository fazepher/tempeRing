
get_HAT_info <- function(mu, l_target, ..., h = .Machine$double.eps^(1/4)){

  mH <- lapply(mu, function(m) -pracma::hessian(l_target, m, h = h, beta = 1, ...))
  Cov <- lapply(mH,solve)
  cholCov <- lapply(Cov,chol)
  detCov <- vapply(cholCov, function(chS) prod(diag(chS))^2, numeric(1))
  l_target_mu <- vapply(mu, l_target, numeric(1), beta = 1, ...)
  w_tilde <- exp(l_target_mu + 0.5*log(detCov))
  w <- w_tilde/sum(w_tilde)

  return(mget(c("mH","Cov","cholCov","detCov","l_target_mu","w","mu")))

}


modAssignment <- function(x, beta, HAT_info){

  lP_vec <-mapply(function(mu,w,Cov) log(w) + mvtnorm::dmvnorm(as.vector(x), mu, Cov/beta, log = TRUE),
                  mu = HAT_info$mu, w = HAT_info$w, Cov = HAT_info$Cov)
  A <- which.max(lP_vec)

  return(list("A" = A, "lP_j" = lP_vec[A]))

}

HAT_rwm_chain <- function(ltemp_target, ..., HAT_info, G_type = "robust",
                          beta_schedule, swap_type = "deo",
                          scale = 1, Temp_Moves = 1000, Within_Moves = 10,
                          x_0 = NULL, x_0_u = NULL, seed = NULL,
                          custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                          silent = FALSE){

  # We create the new target density based on the user-provided target, HAT_info and G_type
  if(G_type == "robust"){

    HAT_l_target <- function(x, beta, HAT_info, ltemp_target, ...){

      ## Basic Weight Preservation

      # Early exit if beta is 1
      if(beta == 1){ return(ltemp_target(x, beta, ...)) }

      # Assign the mode at beta and 1
      mod_beta <- modAssignment(x, beta, HAT_info)
      mod_1 <- modAssignment(x, 1, HAT_info)
      l_mod <- HAT_info$l_target_mu[ mod_beta$A ]

      # Early exit if the assigned modes are the same
      if(mod_beta$A == mod_1$A){ return(ltemp_target(x, beta, ...) + (1-beta) * l_mod) }

      ## G Weight Preservation only when the assigned modes differ

      weight_G <- (2 / (1 + exp(mod_1$lP_j - mod_beta$lP_j))) - 1

      # Early exit if weight_G is negative
      # I DON'T KNOW THE THEORY HERE, PAPER DOESN'T MENTION IT, ASK NICK AND GARETH
      # I WILL RETURN A -Inf log density to always reject a move to such a state
      if(weight_G < 0){
        if(!silent){warning("weight_G is negative")}
        return(-Inf)
      }

      l_cold <- log(1 - weight_G) + ltemp_target(x, 1, ...)
      l_wp_beta <- log(weight_G) + ltemp_target(x, beta, ...) + (1-beta) * l_mod
      l_G <- matrixStats::logSumExp(c(l_cold, l_wp_beta))

      if(!silent && (is.na(l_G) || is.infinite(l_G) || is.nan(l_G))){
        mget(c("l_mod","weight_G","l_cold","l_wp_beta","l_G",
               "x","beta","HAT_info","ltemp_target")) |>
          saveRDS("ALPZ/Error.rds")
      }

      return(l_G)
    }

  }else{

    HAT_l_target <- function(x, beta, HAT_info, ltemp_target, ...){

      ## Basic Weight Preservation

      # Early exit if beta is 1
      if(beta == 1){ return(ltemp_target(x, beta, ...)) }

      # Assign the mode at beta and 1
      mod_beta <- modAssignment(x, beta, HAT_info)
      mod_1 <- modAssignment(x, 1, HAT_info)
      l_mod <- HAT_info$l_target_mu[ mod_beta$A ]

      # Definition depending on whether the assigned modes are the same
      if(mod_beta$A == mod_1$A){
        return(ltemp_target(x, beta, ...) + (1-beta) * l_mod)
      }else{
        # FALTA PREGUNTAR A NICK SI ESTO ESTÁ BIEN (CONFUSA NOTACIÓN EN LOS PAPERS)
        return(0.5*length(x)*(log(2*pi)-log(beta)) + mod_beta$lP_j)
      }

    }

  }


  # Then it's just a wrapper of the PT algorithm for the new HAT target
  PT_rwm_chain(l_target = HAT_l_target, HAT_info = HAT_info, ltemp_target = ltemp_target,
               ...,
               beta_schedule = beta_schedule, swap_type = swap_type,
               scale = scale, Temp_Moves = Temp_Moves, Within_Moves = Within_Moves,
               x_0 = x_0, x_0_u = x_0_u, seed = seed,
               custom_rw_sampler = custom_rw_sampler, target_names = target_names, d = d,
               silent = silent)

}
