
st_temp_step <- function(k_curr, x_curr, l_curr, l_target, ...,
                         beta_schedule, g_schedule = NULL, K = NULL){

  K <- K %||% length(beta_schedule)
  g_schedule <- g_schedule %||% rep(0, K)

  move_right <- runif(1) <= 0.5
  if(k_curr == K && move_right){
    return(list("k_next" = k_curr, "acc" = FALSE, "l_next" = l_curr))
  }
  if(k_curr == 1 && !move_right){
    return(list("k_next" = k_curr, "acc" = FALSE, "l_next" = l_curr))
  }
  k_prop <- ifelse(move_right, k_curr + 1, k_curr - 1)

  delta_g <- g_schedule[k_prop] - g_schedule[k_curr]
  l_prop <- l_target(x_curr, beta_schedule[k_prop], ...)
  delta_H <- l_prop - l_curr + delta_g
  if(delta_H > 0 || log(runif(1)) <= delta_H){
    return(list("k_next" = k_prop, "acc" = TRUE, "l_next" = l_prop))
  }
  return(list("k_next" = k_curr, "acc" = FALSE, "l_next" = l_curr))

}


ST_rwm_chain <- function(l_target, ..., beta_schedule, g_schedule = NULL,
                         scale = 1, Temp_Moves = 1000, Within_Moves = 10,
                         x_0 = NULL, x_0_u = NULL, k_0 = NULL, seed = NULL,
                         custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                         silent = FALSE){

  ## Preparation

  K <- length(beta_schedule)
  g_schedule <- g_schedule %||% rep(0, K)
  # Dimension and Scale
  if(is.list(scale)){
    if(length(scale) != K){
      stop("When scale is a list, it must have the same number of elements as beta_schedule")
    }
    scale_list <- scale
  }else{
    stopifnot(is.numeric(scale))
    scale_list <- rep(list(scale), K)
  }
  d <- d %||% ifelse(is.matrix(scale_list[[1]]), nrow(scale_list[[1]]), length(scale_list[[1]]))
  if(d > 1 && !is.matrix(scale_list[[1]])){
    scale_list <- lapply(scale_list, diag)
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
  }

  # If the user didn't we define proposal sampler(s) as indep. normals
  if(is.null(custom_rw_sampler)){
    sampler <- ifelse(d == 1,
                      function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
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
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }
  k_0 <- k_0 %||% sample(1:K, 1)
  log_lik_0 <- l_target(x_0, beta_schedule[k_0], ...)

  # Preallocate containers
  stopifnot(Temp_Moves > 1)
  stopifnot(Within_Moves >= 1)
  S_Tot <- Temp_Moves*(Within_Moves + 1)
  x <- matrix(nrow = S_Tot, ncol = d, dimnames = list(NULL, target_names))
  acc <- logical(S_Tot)
  log_lik <- numeric(S_Tot)
  k <- numeric(Temp_Moves)


  ## Algorithm

  # First cycle of a Temp_Swap + Within_Moves
  x[1, ] <- x_0
  temp_move <- st_temp_step(k_0, x_0, log_lik_0, l_target, ...,
                            beta_schedule = beta_schedule, g_schedule = g_schedule, K = K)
  log_lik[1] <- temp_move$l_next
  acc[1] <- temp_move$acc
  k[1] <- temp_move$k_next
  kth_sampler <- sampler_list[[ k[1] ]]
  kth_scale <- scale_list[[ k[1] ]]

  rwm_moves <- rwm_sampler_chain(l_target = l_target, beta = beta_schedule[k[1]], ...,
                                 scale = kth_scale, S = Within_Moves,
                                 x_0 = x[1, ], log_lik_0 = log_lik[1],
                                 custom_rw_sampler = kth_sampler)
  window_wm <- 1 + 1:Within_Moves
  x[window_wm, ] <- rwm_moves$x
  log_lik[window_wm] <- rwm_moves$log_lik
  acc[window_wm] <- rwm_moves$acc

  for(c in 2:Temp_Moves){

    s <- c + (c-1)*Within_Moves
    x[s, ] <- x[s-1,]
    temp_move <- st_temp_step(k[c-1], x[s, ], log_lik[s], l_target, ...,
                              beta_schedule = beta_schedule, g_schedule = g_schedule, K = K)
    log_lik[s] <- temp_move$l_next
    acc[s] <- temp_move$acc
    k[c] <- temp_move$k_next
    kth_sampler <- sampler_list[[ k[c] ]]
    kth_scale <- scale_list[[ k[c] ]]

    rwm_moves <- rwm_sampler_chain(l_target = l_target, beta = beta_schedule[k[c]], ...,
                                   scale = kth_scale, S = Within_Moves,
                                   x_0 = x[s, ], log_lik_0 = log_lik[s],
                                   custom_rw_sampler = kth_sampler)
    window_wm <- s + 1:(Within_Moves)
    x[window_wm, ] <- rwm_moves$x
    log_lik[window_wm] <- rwm_moves$log_lik
    acc[window_wm] <- rwm_moves$acc

  }

  return(mget(c("x","k","acc","log_lik")))

}
