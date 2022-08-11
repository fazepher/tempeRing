
PT_rwm_chain <- function(l_target, ..., beta_schedule, swap_type = "deo",
                         scale = 1, Temp_Moves = 1000, Within_Moves = 10,
                         x_0 = NULL, x_0_u = NULL, seed = NULL,
                         custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                         silent = FALSE){


  stopifnot(swap_type %in% c("deo","seo","naive"))
  K <- length(beta_schedule)
  if(swap_type != "naive"){
    odd_indices <- seq(1, K, by = 2)
    even_indices <- seq(2, K, by = 2)
  }

  # Dimension and Scales
  if(is.list(scale)){
    if(length(scale) != K){
      stop("As a list, length of scale must must match that of beta_schedule")
    }
    scale_list <- scale
  }else{
    stopifnot(is.numeric(scale))
    scale_list <- lapply(beta_schedule, function(beta) scale/sqrt(beta))
  }
  d <- d %||% ifelse(is.matrix(scale_list[[1]]), nrow(scale_list[[1]]), length(scale_list[[1]]))
  if(d > 1 && !is.matrix(scale_list[[1]])){
    scale_list <- lapply(scale_list, function(scale) diag(scale, d))
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
  }

  # If the user didn't we define proposal sampler(s) as indep. normals
  if(is.null(custom_rw_sampler)){
    sampler <- ifelse(d == 1,
                      function(x, scale){ stats::rnorm(n = 1, mean = x, sd = scale) },
                      function(x, scale){ mvtnorm::rmvnorm(n = 1, mean = x, sigma = scale) })
    sampler_list <- rep(list(sampler), K)
  }else{
    sampler_list <- ifelse(is.list(custom_rw_sampler),
                           custom_rw_sampler,
                           rep(list(custom_rw_sampler),K))
  }

  # Starting point
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(x_0)){
    x_0_u <- x_0_u %||% 2
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- matrix(data = stats::runif(d*K, min = -x_0_u, max = x_0_u), nrow = K, ncol = d)
  }else{
    stopifnot(is.numeric(x_0))
    stopifnot(nrow(x_0) == K && ncol(x_0) == d)
  }
  log_lik_0 <- mapply(function(x,beta) l_target(x, beta, ...),
                      asplit(x_0, 1),  beta_schedule)

  # Preallocate containers
  stopifnot(Temp_Moves > 1)
  stopifnot(Within_Moves >= 1)
  S_Tot <- Temp_Moves*(Within_Moves + 1)
  x <- array(dim = c(S_Tot, K, d), dimnames = list(NULL,NULL,target_names))
  acc <- matrix(nrow = S_Tot, ncol = K)
  log_lik <- matrix(nrow = S_Tot, ncol = K)
  k_indexes <- matrix(nrow = Temp_Moves, ncol = K)
  beta_indexes <- matrix(nrow = Temp_Moves, ncol = K)

  ## Algorithm
  window_wm <- 1 + 1:(Within_Moves)
  # First Cycle
  swap_move <- temp_swap_move(type = swap_type, c = 1,
                              array(x_0, dim = c(1, K, d)),
                              beta_schedule,
                              1:K,
                              log_lik_0,
                              l_target, ...,
                              K = K, d = d)
  x[1, , ] <- x_0
  log_lik[1, ] <- swap_move$l_next
  acc[1, ] <- swap_move$acc
  beta_indexes[1, ] <- swap_move$beta_next
  k_indexes[1, ] <- swap_move$k_next

  for(k in 1:K){
    rwm_moves <- rwm_sampler_chain(l_target = l_target,
                                   beta = beta_indexes[1, k],
                                   ...,
                                   scale = scale_list[[ k_indexes[1, k] ]],
                                   S = Within_Moves,
                                   x_0 = x[1, k, ],
                                   log_lik_0 = log_lik[1, k],
                                   custom_rw_sampler = sampler_list[[ k_indexes[1, k] ]])
    x[window_wm, k, ] <- rwm_moves$x
    acc[window_wm, k] <- rwm_moves$acc
    log_lik[window_wm, k] <- rwm_moves$log_lik
  }

  # Next Cycles
  for(c in 2:Temp_Moves){

    s <- c + (c-1)*Within_Moves
    window_wm <- s + 1:(Within_Moves)
    swap_move <- temp_swap_move(type = swap_type, c = c,
                                x[s-1, , , drop = FALSE],
                                beta_indexes[c-1, ],
                                k_indexes[c-1, ],
                                log_lik[s-1, ],
                                l_target, ...,
                                K = K, d = d)
    x[s, , ] <- x[s-1, , , drop = FALSE]
    log_lik[s, ] <- swap_move$l_next
    acc[s, ] <- swap_move$acc
    beta_indexes[c, ] <- swap_move$beta_next
    k_indexes[c, ] <- swap_move$k_next

    for(k in 1:K){
      rwm_moves <- rwm_sampler_chain(l_target = l_target,
                                     beta = beta_indexes[c, k],
                                     ...,
                                     scale = scale_list[[ k_indexes[c, k] ]],
                                     S = Within_Moves,
                                     x_0 = x[s, k, ],
                                     log_lik_0 = log_lik[s, k],
                                     custom_rw_sampler = sampler_list[[ k_indexes[c, k] ]])
      x[window_wm, k, ] <- rwm_moves$x
      acc[window_wm, k] <- rwm_moves$acc
      log_lik[window_wm, k] <- rwm_moves$log_lik
    }

  }

  return(mget(c("x","acc","log_lik","k_indexes","beta_indexes")))

}
