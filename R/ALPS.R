
#' @export
ALPS_rwm_chain <- function(ltemp_target, ..., HAT = TRUE, HAT_info,
                           beta_schedule, swap_type = "deo", quanta_levels = NULL,
                           scale = 1, Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                           x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, jump_p = 0.7,
                           k_0 = NULL, custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                           quanta_mode_info = NULL, silent = FALSE){

  #--- HAT use -------------
  l_target <- lHAT_target
  target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                   rlang::dots_list(...))


  #--- Preparation -------------

  # General PT parameters
  start_time <- Sys.time()
  global_times <- rep(start_time, 3)
  K <- length(beta_schedule)
  odd_indices <- seq(1, K, by = 2)
  even_indices <- seq(2, K, by = 2)

  # Leap-Point Sampler
  lpsampler <- function(x_curr, beta_max = beta_schedule[K], mode_info = HAT_info){
    z <- sample(x = seq_along(mode_info$w),size = 1, prob = mode_info$w)
    x_prop <- rmvtnorm_temp(n = 1, beta = beta_max, mu = mode_info$modes[[z]],
                            LChol_sigma = t(mode_info$cholCov[[z]]))
    return(x_prop)
  }
  lpsampler_q <- function(x, beta_max = beta_schedule[K], mode_info = HAT_info){
    n_modes <- length(mode_info$w)
    l_modes <- vapply(1:n_modes, function(m)
      lmvtnorm_temp(x = x, beta = beta_max, mu = mode_info$modes[[m]],
                    sigma_inv = mode_info$mH[[m]], logdet_sigma = 2*mode_info$half_l_detCov[[m]]) +
        log(mode_info$w[[m]]), FUN.VALUE = 1.0)

    indmax <- which.max(l_modes)   ##### Create a stable evaluation of the log density.
    comb <- exp(l_modes[-indmax]-l_modes[indmax])
    return( l_modes[indmax] + log(1+sum(comb)) )

  }


  # QuanTA levels checking
  if(!is.null(quanta_levels)){
    if(!all(quanta_levels %in% seq(1,K-1))){
      stop(paste("Invalid quanta_levels.",
                 "These must be the base indexes for which QuanTA swaps will be attempted,",
                 "hence a vector of integers between 1 and K-1,",
                 "where K is the number of temperatures used.",
                 sep = "\n"))
    }
    quanta_mode_info <- quanta_mode_info %||% HAT_info
  }

  # Dimension and Scales
  if(is.list(scale)){
    if(length(scale) != K){
      stop("As a list, length of scale must must equal that of beta_schedule")
    }
    scale_list <- scale
  }else{
    stopifnot(is.numeric(scale))
    scale_list <- lapply(beta_schedule, function(beta) scale/beta)
  }
  d <- d %||% ifelse(is.matrix(scale_list[[1]]), nrow(scale_list[[1]]), length(scale_list[[1]]))
  if(d > 1 && !is.matrix(scale_list[[1]])){
    scale_list <- lapply(scale_list, function(scale) diag(scale, d))
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
  }

  # Checking for valid sample sizes
  stopifnot(Cycles >= 1)
  stopifnot(-1 <= burn_cycles && burn_cycles < Cycles)
  stopifnot(Temp_Moves >= 1)
  stopifnot(Within_Moves >= 1)

  # If the user didn't we define proposal sampler(s) as indep. normals
  if(is.list(custom_rw_sampler)){
    if(length(custom_rw_sampler) != K){
      stop("As a list, custom_rw_sampler must have the same length of beta_schedule")
    }
    sampler_list <- custom_rw_sampler
  }else{
    sampler <- custom_rw_sampler %||%
      ifelse(d == 1,
             function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
             function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })
    sampler_list <- rep(list(sampler),K)
  }

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- matrix(data = runif(d*K, min = -x_0_u, max = x_0_u), nrow = K, ncol = d)
  }else{
    stopifnot(is.numeric(x_0))
    stopifnot(nrow(x_0) == K && ncol(x_0) == d)
  }
  k_0 <- k_0 %||% 1:K
  b_0 <- beta_schedule[k_0]
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(m in 1:K){
      l_0[m] <- do.call(l_target, c(list(x = x_0[k_0[m], ], beta = b_0[m]), target_args))
    }
  }

  # Preallocate containers
  cycle_length <- Temp_Moves + Within_Moves
  S_Tot <- Cycles*(cycle_length)
  x <- array(dim = c(S_Tot + 1, K, d), dimnames = list(NULL, NULL, target_names))
  l_x <- matrix(nrow = S_Tot + 1, ncol = K)
  y <- array(dim = c(S_Tot, K, d), dimnames = list(NULL, NULL, target_names))
  l_y <- matrix(nrow = S_Tot, ncol = K)
  delta_l <- matrix(nrow = S_Tot, ncol = K)
  k_indexes <- array(dim = c(Cycles + 1, Temp_Moves + 1, K))
  beta_indexes <- array(dim = c(Cycles + 1, Temp_Moves + 1, K))
  swap_acc <- array(-1, dim = c(Cycles, Temp_Moves, K))
  rwm_acc <- array(dim = c(Cycles, Within_Moves, K))
  cycle_times <- as.POSIXct(rep(NA, Cycles))

  #--- Algorithm -------------

  # Initialize
  k_indexes[1, 1, ] <- k_0
  beta_indexes[1, 1, ] <- b_0[k_0]
  x[1, , ] <- x_0
  l_x[1, ] <- l_0

  # Determine LPS cycles at coldest level
  u_lps <- runif(Cycles) <= jump_p

  # Run iterations
  for(c in 1:Cycles){

    cycle_times[c] <- Sys.time()
    if(!silent & isTRUE(c %% floor(Cycles*0.05) == 0)){
      cat(paste0("Avance: ",round(100*c/Cycles),"%"),sep = "\n")
      print(Sys.time())
    }

    # Current cycle position index
    i <- (c-1)*(cycle_length) + 2

    # Temperature Swaps
    for(j in 1:Temp_Moves){

      swap_arguments <- list(type = swap_type, j_deo = ifelse(Temp_Moves == 1, c, j),
                             quanta_levels = quanta_levels, mode_info = quanta_mode_info,
                             beta_curr = beta_indexes[c, j, ],
                             k_curr = k_indexes[c, j,  ],
                             x_curr = x[(i-2)+j, , , drop = FALSE],
                             l_curr = l_x[(i-2)+j, ],
                             l_target, target_args,
                             K = K,
                             odd_indices = odd_indices,
                             even_indices = even_indices,
                             d = d)
      swap_move <- do.call(alps_swap_move, swap_arguments)

      swap_acc[c, j, ] <- swap_move$acc
      k_indexes[c, j+1, ] <- swap_move$k_next
      beta_indexes[c, j+1, ] <- swap_move$beta_next

      next_ind <- (i-1) + j
      prop_ind <- next_ind - 1

      x[next_ind, , ] <- swap_move$x_next
      l_x[next_ind, ] <- swap_move$l_next
      y[prop_ind, , ] <- swap_move$x_prop
      l_y[prop_ind, ] <- swap_move$l_prop

    }

    # Within temperature moves

    # Update current cycle position index
    i <- i + Temp_Moves - 1
    within_window_next <- i + 1:Within_Moves
    within_window_prop <- within_window_next - 1

    # Random Walk Metropolis
    for(k in 1:(K-1)){

      k_i <- which(k_indexes[c, Temp_Moves + 1, ] == k)

      rwm_level_args <- list(x_0 = x[i, k_i , ],
                             l_0 = l_x[i, k_i],
                             beta = beta_schedule[k],
                             custom_rw_sampler = sampler_list[[k]],
                             scale = scale_list[[k]],
                             S = Within_Moves, burn = 0, silent = TRUE)
      rwm_moves <- do.call(rwm_sampler_chain, c(rwm_level_args, l_target, target_args))

      rwm_acc[c, , k] <- rwm_moves$acc
      x[within_window_next, k_i, ] <- rwm_moves$x
      l_x[within_window_next, k_i] <- rwm_moves$l_x
      y[within_window_prop, k_i, ] <- rwm_moves$y
      l_y[within_window_prop, k_i] <- rwm_moves$l_y
      delta_l[within_window_prop, k_i] <- rwm_moves$delta_l

    }
    # Leap Sampler at Coldest Level
    k_i <- which(k_indexes[c, Temp_Moves + 1, ] == K)
    if(u_lps[c]){

      for(s in 1:Within_Moves){
        x_curr_lps <- x[i+s-1, k_i , ]
        l_curr_lps <- l_x[i+s-1, k_i]
        x_prop_lps <- lpsampler(x_curr = x_curr_lps)
        l_prop_lps <- do.call(l_target, c(list(x = x_prop_lps, beta = beta_schedule[K]), target_args))
        lsaq_c2p <- lpsampler_q(x = x_prop_lps)
        lsaq_p2c <- lpsampler_q(x = x_curr_lps)
        lsa_step <- mh_step(x_curr = x_curr_lps, x_prop = x_prop_lps,
                            l_curr = l_curr_lps, l_prop = l_prop_lps,
                            lq_c2p = lsaq_c2p, lq_p2c = lsaq_p2c,
                            do_checks = FALSE)

        rwm_acc[c, s, K] <- lsa_step$accepted
        ww_n <- i + s
        ww_p <- ww_n - 1
        x[ww_n, k_i, ] <- lsa_step$x_next
        l_x[ww_n, k_i] <- lsa_step$l_next
        y[ww_p, k_i, ] <- lsa_step$x_prop
        l_y[ww_p, k_i] <- lsa_step$l_prop
        delta_l[ww_p, k_i] <- lsa_step$delta_l
      }

    }else{

      rwm_level_args <- list(x_0 = x[i, k_i , ],
                             l_0 = l_x[i, k_i],
                             beta = beta_schedule[K],
                             custom_rw_sampler = sampler_list[[K]],
                             scale = scale_list[[K]],
                             S = Within_Moves, burn = 0, silent = TRUE)
      rwm_moves <- do.call(rwm_sampler_chain, c(rwm_level_args, l_target, target_args))

      rwm_acc[c, , K] <- rwm_moves$acc
      x[within_window_next, k_i, ] <- rwm_moves$x
      l_x[within_window_next, k_i] <- rwm_moves$l_x
      y[within_window_prop, k_i, ] <- rwm_moves$y
      l_y[within_window_prop, k_i] <- rwm_moves$l_y
      delta_l[within_window_prop, k_i] <- rwm_moves$delta_l

    }


    # Copy parameters for next cycle
    k_indexes[c+1, 1, ] <- k_indexes[c, Temp_Moves + 1, ]
    beta_indexes[c+1, 1, ] <- beta_indexes[c, Temp_Moves + 1, ]

  }
  global_times[2] <- Sys.time()

  #--- Post processing and result -------------

  swap_acc_rates <- apply(swap_acc, 3, function(swaps) mean(swaps[swaps>=0]))
  rwm_acc_rates <- colMeans(rwm_acc, dims = 2)

  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Swap Acceptance Rates:", round(swap_acc_rates,3)), sep = "\n")
    cat(paste("RWM Acceptance Rates:", round(rwm_acc_rates,3)), sep = "\n")
  }

  if(burn_cycles == -1){
    global_times[3] <- Sys.time()
    return(list(x = x, l_x = l_x, y = y, l_y = l_y, delta_l = delta_l, u_lps = u_lps,
                k_indexes = k_indexes[-(Cycles+1), , ],
                beta_indexes = beta_indexes[-(Cycles+1), , ],
                swap_acc = swap_acc, swap_acc_rates = swap_acc_rates,
                rwm_acc = rwm_acc, rwm_acc_rates = rwm_acc_rates,
                cycle_times = cycle_times, global_times = global_times))
  }

  burn_window <- seq(1,burn_cycles*cycle_length + 1)
  x_r <- x[-burn_window, , , drop = FALSE]
  l_x_r <- l_x[-burn_window, ]
  y_r <- y[-burn_window, , , drop = FALSE]
  l_y_r <- l_y[-burn_window, ]
  delta_l_r <- delta_l[-burn_window, ]
  u_lps_r <- u_lps[-burn_window, ]

  if(burn_cycles == 0){
    rwm_acc_r <- rwm_acc
    swap_acc_r <- swap_acc
    b_r <- beta_indexes[-c(1,Cycles+1), , ]
    k_r <- k_indexes[-c(1,Cycles+1), , ]
  }else{
    burn_window <- 1:burn_cycles
    rwm_acc_r <- rwm_acc[-burn_window, , ]
    swap_acc_r <- swap_acc[-burn_window, ]
    b_r <- beta_indexes[-c(burn_window,Cycles + 1), , ]
    k_r <- k_indexes[-c(burn_window,Cycles + 1), , ]
  }

  global_times[3] <- Sys.time()
  return(list(x = x_r, l_x = l_x_r, y = y_r, l_y = l_y_r, delta_l = delta_l_r, u_lps = u_lps_r,
              k_indexes = k_r, beta_indexes = b_r,
              swap_acc = swap_acc_r, swap_acc_rates = swap_acc_rates,
              rwm_acc = rwm_acc_r, rwm_acc_rates = rwm_acc_rates,
              cycle_times = cycle_times, global_times = global_times))

}


#' @export
ALPS_rwm_leaner_chain <- function(ltemp_target, ..., HAT = TRUE, HAT_info,
                                  beta_schedule, swap_type = "deo", quanta_levels = NULL,
                                  scale = 1, Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                                  x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, jump_p = 0.7,
                                  k_0 = NULL, custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                                  quanta_mode_info = NULL, silent = FALSE){

  #--- HAT use -------------
  l_target <- lHAT_target
  target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                   rlang::dots_list(...))

  #--- Preparation -------------

  # General PT parameters
  start_time <- Sys.time()
  global_times <- rep(start_time, 3)
  K <- length(beta_schedule)
  odd_indices <- seq(1, K, by = 2)
  even_indices <- seq(2, K, by = 2)

  # Leap-Point Sampler
  lpsampler <- function(x_curr, beta_max = beta_schedule[K], mode_info = HAT_info){
    z <- sample(x = seq_along(mode_info$w),size = 1, prob = mode_info$w)
    x_prop <- rmvtnorm_temp(n = 1, beta = beta_max, mu = mode_info$modes[[z]],
                            LChol_sigma = t(mode_info$cholCov[[z]]))
    return(x_prop)
  }
  lpsampler_q <- function(x, beta_max = beta_schedule[K], mode_info = HAT_info){
    n_modes <- length(mode_info$w)
    l_modes <- vapply(1:n_modes, function(m)
      lmvtnorm_temp(x = x, beta = beta_max, mu = mode_info$modes[[m]],
                    sigma_inv = mode_info$mH[[m]], logdet_sigma = 2*mode_info$half_l_detCov[[m]]) +
        log(mode_info$w[[m]]), FUN.VALUE = 1.0)

    indmax <- which.max(l_modes)   ##### Create a stable evaluation of the log density.
    comb <- exp(l_modes[-indmax]-l_modes[indmax])
    return( l_modes[indmax] + log(1+sum(comb)) )

  }


  # QuanTA levels checking
  if(is.null(quanta_levels)){
    quanta_levels <- K-1
  }else{
    if(!all(quanta_levels %in% seq(1,K-1))){
      stop(paste("Invalid quanta_levels.",
                 "These must be the base indexes for which QuanTA swaps will be attempted,",
                 "hence a vector of integers between 1 and K-1,",
                 "where K is the number of temperatures used.",
                 sep = "\n"))
    }
    quanta_mode_info <- quanta_mode_info %||% HAT_info
  }

  # Dimension and Scales
  if(is.list(scale)){
    if(length(scale) != K){
      stop("As a list, length of scale must must equal that of beta_schedule")
    }
    scale_list <- scale
  }else{
    stopifnot(is.numeric(scale))
    scale_list <- lapply(beta_schedule, function(beta) scale/beta)
  }
  d <- d %||% ifelse(is.matrix(scale_list[[1]]), nrow(scale_list[[1]]), length(scale_list[[1]]))
  if(d > 1 && !is.matrix(scale_list[[1]])){
    scale_list <- lapply(scale_list, function(scale) diag(scale, d))
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
  }

  # Checking for valid sample sizes
  stopifnot(Cycles >= 1)
  stopifnot(-1 <= burn_cycles && burn_cycles < Cycles)
  stopifnot(Temp_Moves >= 1)
  stopifnot(Within_Moves >= 1)

  # If the user didn't we define proposal sampler(s) as indep. normals
  if(is.list(custom_rw_sampler)){
    if(length(custom_rw_sampler) != K){
      stop("As a list, custom_rw_sampler must have the same length of beta_schedule")
    }
    sampler_list <- custom_rw_sampler
  }else{
    sampler <- custom_rw_sampler %||%
      ifelse(d == 1,
             function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
             function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })
    sampler_list <- rep(list(sampler),K)
  }

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- matrix(data = runif(d*K, min = -x_0_u, max = x_0_u), nrow = K, ncol = d)
  }else{
    stopifnot(is.numeric(x_0))
    stopifnot(nrow(x_0) == K && ncol(x_0) == d)
  }
  k_0 <- k_0 %||% 1:K
  b_0 <- beta_schedule[k_0]
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(m in 1:K){
      l_0[m] <- do.call(l_target, c(list(x = x_0[k_0[m], ], beta = b_0[m]), target_args))
    }
  }

  # Preallocate containers
  cycle_length <- Temp_Moves + Within_Moves
  S_Tot <- Cycles*(cycle_length)
  x <- array(dim = c(S_Tot + 1, K, d), dimnames = list(NULL, NULL, target_names))
  k_indexes <- array(dim = c(Cycles + 1, Temp_Moves + 1, K))
  beta_indexes <- array(dim = c(Cycles + 1, Temp_Moves + 1, K))
  swap_acc <- array(-1, dim = c(Cycles, Temp_Moves, K))
  rwm_acc <- array(dim = c(Cycles, Within_Moves, K))
  cycle_times <- as.POSIXct(rep(NA, Cycles))

  #--- Algorithm -------------

  # Initialize
  k_indexes[1, 1, ] <- k_0
  beta_indexes[1, 1, ] <- b_0[k_0]
  x[1, , ] <- x_0
  l_x <- l_0

  # Determine LPS cycles at coldest level
  u_lps <- runif(Cycles) <= jump_p

  # Run iterations
  for(c in 1:Cycles){

    cycle_times[c] <- Sys.time()
    if(!silent & isTRUE(c %% floor(Cycles*0.05) == 0)){
      cat(paste0("Avance: ",round(100*c/Cycles),"%"),sep = "\n")
      print(Sys.time())
    }

    # Current cycle position index
    i <- (c-1)*(cycle_length) + 2

    # Temperature Swaps
    for(j in 1:Temp_Moves){

      swap_arguments <- list(type = swap_type, j_deo = ifelse(Temp_Moves == 1, c, j),
                             quanta_levels = quanta_levels, mode_info = quanta_mode_info,
                             beta_curr = beta_indexes[c, j, ],
                             k_curr = k_indexes[c, j,  ],
                             x_curr = x[(i-2)+j, , , drop = FALSE],
                             l_curr = l_x,
                             l_target, target_args,
                             K = K,
                             odd_indices = odd_indices,
                             even_indices = even_indices,
                             d = d)
      swap_move <- do.call(alps_swap_move, swap_arguments)

      swap_acc[c, j, ] <- swap_move$acc
      k_indexes[c, j+1, ] <- swap_move$k_next
      beta_indexes[c, j+1, ] <- swap_move$beta_next

      next_ind <- (i-1) + j
      prop_ind <- next_ind - 1

      x[next_ind, , ] <- swap_move$x_next
      l_x <- swap_move$l_next

    }

    # Within temperature moves

    # Update current cycle position index
    i <- i + Temp_Moves - 1
    within_window_next <- i + 1:Within_Moves
    within_window_prop <- within_window_next - 1

    # Random Walk Metropolis
    for(k in 1:(K-1)){

      k_i <- which(k_indexes[c, Temp_Moves + 1, ] == k)

      rwm_level_args <- list(x_0 = x[i, k_i , ],
                             l_0 = l_x[k_i],
                             beta = beta_schedule[k],
                             custom_rw_sampler = sampler_list[[k]],
                             scale = scale_list[[k]],
                             S = Within_Moves, burn = 0, silent = TRUE)
      rwm_moves <- do.call(rwm_sampler_leaner_chain, c(rwm_level_args, l_target, target_args))

      rwm_acc[c, , k] <- rwm_moves$acc
      x[within_window_next, k_i, ] <- rwm_moves$x
      l_x[k_i] <- rwm_moves$l_x_curr

    }
    # Leap Sampler at Coldest Level
    k_i <- which(k_indexes[c, Temp_Moves + 1, ] == K)
    if(u_lps[c]){

      for(s in 1:Within_Moves){
        x_curr_lps <- x[i+s-1, k_i , ]
        l_curr_lps <- l_x[k_i]
        x_prop_lps <- lpsampler(x_curr = x_curr_lps)
        l_prop_lps <- do.call(l_target, c(list(x = x_prop_lps, beta = beta_schedule[K]), target_args))
        lsaq_c2p <- lpsampler_q(x = x_prop_lps)
        lsaq_p2c <- lpsampler_q(x = x_curr_lps)
        lsa_step <- mh_step(x_curr = x_curr_lps, x_prop = x_prop_lps,
                            l_curr = l_curr_lps, l_prop = l_prop_lps,
                            lq_c2p = lsaq_c2p, lq_p2c = lsaq_p2c,
                            do_checks = FALSE, full_return = FALSE)

        rwm_acc[c, s, K] <- lsa_step$accepted
        ww_n <- i + s
        ww_p <- ww_n - 1
        x[ww_n, k_i, ] <- lsa_step$x_next
        l_x[k_i] <- lsa_step$l_next
      }

    }else{

      rwm_level_args <- list(x_0 = x[i, k_i , ],
                             l_0 = l_x[k_i],
                             beta = beta_schedule[K],
                             custom_rw_sampler = sampler_list[[K]],
                             scale = scale_list[[K]],
                             S = Within_Moves, burn = 0, silent = TRUE)
      rwm_moves <- do.call(rwm_sampler_chain, c(rwm_level_args, l_target, target_args))

      rwm_acc[c, , K] <- rwm_moves$acc
      x[within_window_next, k_i, ] <- rwm_moves$x
      l_x[k_i] <- rwm_moves$l_x[Within_Moves]

    }


    # Copy parameters for next cycle
    k_indexes[c+1, 1, ] <- k_indexes[c, Temp_Moves + 1, ]
    beta_indexes[c+1, 1, ] <- beta_indexes[c, Temp_Moves + 1, ]

  }
  global_times[2] <- Sys.time()

  #--- Post processing and result -------------

  swap_acc_rates <- apply(swap_acc, 3, function(swaps) mean(swaps[swaps>=0]))
  rwm_acc_rates <- colMeans(rwm_acc, dims = 2)

  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Swap Acceptance Rates:", round(swap_acc_rates,3)), sep = "\n")
    cat(paste("RWM Acceptance Rates:", round(rwm_acc_rates,3)), sep = "\n")
  }

  if(burn_cycles == -1){
    global_times[3] <- Sys.time()
    return(list(x = x,
                k_indexes = k_indexes[-(Cycles+1), , ],
                beta_indexes = beta_indexes[-(Cycles+1), , ],
                swap_acc = swap_acc, swap_acc_rates = swap_acc_rates,
                rwm_acc = rwm_acc, rwm_acc_rates = rwm_acc_rates,
                cycle_times = cycle_times, global_times = global_times))
  }

  burn_window <- seq(1,burn_cycles*cycle_length + 1)
  x_r <- x[-burn_window, , , drop = FALSE]

  if(burn_cycles == 0){
    rwm_acc_r <- rwm_acc
    swap_acc_r <- swap_acc
    b_r <- beta_indexes[-c(1,Cycles+1), , ]
    k_r <- k_indexes[-c(1,Cycles+1), , ]
  }else{
    burn_window <- 1:burn_cycles
    rwm_acc_r <- rwm_acc[-burn_window, , ]
    swap_acc_r <- swap_acc[-burn_window, ]
    b_r <- beta_indexes[-c(burn_window,Cycles + 1), , ]
    k_r <- k_indexes[-c(burn_window,Cycles + 1), , ]
  }

  global_times[3] <- Sys.time()
  return(list(x = x_r,
              k_indexes = k_r, beta_indexes = b_r,
              swap_acc = swap_acc_r, swap_acc_rates = swap_acc_rates,
              rwm_acc = rwm_acc_r, rwm_acc_rates = rwm_acc_rates,
              cycle_times = cycle_times, global_times = global_times))

}

#' @export
ALPS_mh_chain <- function(ltemp_target, ..., d,
                          sampler_list, sampler_args_list = NULL,
                          lq_list, lq_args_list = NULL,
                          HAT = TRUE, HAT_info,
                          beta_schedule, swap_type = "deo", quanta_levels = NULL,
                          Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                          x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, jump_p = 0.7,
                          k_0 = NULL, quanta_mode_info = NULL, silent = FALSE, target_names = NULL){

  #--- HAT use -------------
  l_target <- lHAT_target
  target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                   rlang::dots_list(...))

  #--- Preparation -------------

  # General PT parameters
  start_time <- Sys.time()
  global_times <- rep(start_time, 3)
  K <- length(beta_schedule)
  odd_indices <- seq(1, K, by = 2)
  even_indices <- seq(2, K, by = 2)

  # Leap-Point Sampler
  lpsampler <- function(x_curr, beta_max = beta_schedule[K], mode_info = HAT_info){
    z <- sample(x = seq_along(mode_info$w),size = 1, prob = mode_info$w)
    x_prop <- rmvtnorm_temp(n = 1, beta = beta_max, mu = mode_info$modes[[z]],
                            LChol_sigma = t(mode_info$cholCov[[z]]))
    return(x_prop)
  }
  lpsampler_q <- function(x, beta_max = beta_schedule[K], mode_info = HAT_info){
    n_modes <- length(mode_info$w)
    l_modes <- vapply(1:n_modes, function(m)
      lmvtnorm_temp(x = x, beta = beta_max, mu = mode_info$modes[[m]],
                    sigma_inv = mode_info$mH[[m]], logdet_sigma = 2*mode_info$half_l_detCov[[m]]) +
        log(mode_info$w[[m]]), FUN.VALUE = 1.0)

    indmax <- which.max(l_modes)   ##### Create a stable evaluation of the log density.
    comb <- exp(l_modes[-indmax]-l_modes[indmax])
    return( l_modes[indmax] + log(1+sum(comb)) )

  }


  # QuanTA levels checking
  if(is.null(quanta_levels)){
    quanta_levels <- K-1
  }else{
    if(!all(quanta_levels %in% seq(1,K-1))){
      stop(paste("Invalid quanta_levels.",
                 "These must be the base indexes for which QuanTA swaps will be attempted,",
                 "hence a vector of integers between 1 and K-1,",
                 "where K is the number of temperatures used.",
                 sep = "\n"))
    }
    quanta_mode_info <- quanta_mode_info %||% HAT_info
  }


  # Checking for valid sample sizes
  stopifnot(Cycles >= 1)
  stopifnot(-1 <= burn_cycles && burn_cycles < Cycles)
  stopifnot(Temp_Moves >= 1)
  stopifnot(Within_Moves >= 1)

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- matrix(data = runif(d*K, min = -x_0_u, max = x_0_u), nrow = K, ncol = d)
  }else{
    stopifnot(is.numeric(x_0))
    stopifnot(nrow(x_0) == K && ncol(x_0) == d)
  }
  k_0 <- k_0 %||% 1:K
  b_0 <- beta_schedule[k_0]
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(m in 1:K){
      l_0[m] <- do.call(l_target, c(list(x = x_0[k_0[m], ], beta = b_0[m]), target_args))
    }
  }

  # Preallocate containers
  cycle_length <- Temp_Moves + Within_Moves
  S_Tot <- Cycles*(cycle_length)
  x <- array(dim = c(S_Tot + 1, K, d), dimnames = list(NULL, NULL, target_names))
  l_x <- matrix(nrow = S_Tot + 1, ncol = K)
  y <- array(dim = c(S_Tot, K, d), dimnames = list(NULL, NULL, target_names))
  l_y <- matrix(nrow = S_Tot, ncol = K)
  delta_l <- matrix(nrow = S_Tot, ncol = K)
  k_indexes <- array(dim = c(Cycles + 1, Temp_Moves + 1, K))
  beta_indexes <- array(dim = c(Cycles + 1, Temp_Moves + 1, K))
  swap_acc <- array(-1, dim = c(Cycles, Temp_Moves, K))
  mh_acc <- array(dim = c(Cycles, Within_Moves, K))
  cycle_times <- as.POSIXct(rep(NA, Cycles))

#--- Algorithm -------------

  # Initialize
  k_indexes[1, 1, ] <- k_0
  beta_indexes[1, 1, ] <- b_0
  x[1, , ] <- x_0
  l_x[1, ] <- l_0

  # Determine LPS cycles at coldest level
  u_lps <- runif(Cycles) <= jump_p

  # Run iterations
  for(c in 1:Cycles){

    cycle_times[c] <- Sys.time()
    if(!silent & isTRUE(c %% floor(Cycles*0.05) == 0)){
      cat(paste0("Avance: ",round(100*c/Cycles),"%"),sep = "\n")
      print(Sys.time())
    }

    # Current cycle position index
    i <- (c-1)*(cycle_length) + 2

    # Temperature Swaps
    for(j in 1:Temp_Moves){

      swap_arguments <- list(type = swap_type, j_deo = ifelse(Temp_Moves == 1, c, j),
                             quanta_levels = quanta_levels, mode_info = quanta_mode_info,
                             beta_curr = beta_indexes[c, j, ],
                             k_curr = k_indexes[c, j,  ],
                             x_curr = x[(i-2)+j, , , drop = FALSE],
                             l_curr = l_x[(i-2)+j, ],
                             l_target, target_args,
                             K = K,
                             odd_indices = odd_indices,
                             even_indices = even_indices,
                             d = d)
      swap_move <- do.call(alps_swap_move, swap_arguments)

      swap_acc[c, j, ] <- swap_move$acc
      k_indexes[c, j+1, ] <- swap_move$k_next
      beta_indexes[c, j+1, ] <- swap_move$beta_next

      next_ind <- (i-1) + j
      prop_ind <- next_ind - 1

      x[next_ind, , ] <- swap_move$x_next
      l_x[next_ind, ] <- swap_move$l_next
      y[prop_ind, , ] <- swap_move$x_prop
      l_y[prop_ind, ] <- swap_move$l_prop

    }

    # Within temperature moves

    # Update current cycle position index
    i <- i + Temp_Moves - 1
    within_window_next <- i + 1:Within_Moves
    within_window_prop <- within_window_next - 1

    # Random Walk Metropolis
    for(k in 1:(K-1)){

      k_i <- which(k_indexes[c, Temp_Moves + 1, ] == k)

      mh_level_args <- list(x_0 = x[i, k_i , ],
                            l_0 = l_x[i, k_i],
                            beta = beta_schedule[k],
                            mh_sampler = sampler_list[[k]],
                            other_sampler_args = sampler_args_list[[k]],
                            lq_mh = lq_list[[k]],
                            other_lq_args = lq_args_list[[k]],
                            S = Within_Moves, burn = 0, silent = TRUE, d = d)
      mh_moves <- do.call(mh_sampler_chain, c(mh_level_args,l_target, target_args))

      mh_acc[c, , k] <- mh_moves$acc
      x[within_window_next, k_i, ] <- mh_moves$x
      l_x[within_window_next, k_i] <- mh_moves$l_x
      y[within_window_prop, k_i, ] <- mh_moves$y
      l_y[within_window_prop, k_i] <- mh_moves$l_y
      delta_l[within_window_prop, k_i] <- mh_moves$delta_l

    }
    # Leap Sampler at Coldest Level
    k_i <- which(k_indexes[c, Temp_Moves + 1, ] == K)
    if(u_lps[c]){

      for(s in 1:Within_Moves){
        x_curr_lps <- x[i+s-1, k_i , ]
        l_curr_lps <- l_x[i+s-1, k_i]
        x_prop_lps <- lpsampler(x_curr = x_curr_lps)
        l_prop_lps <- do.call(l_target, c(list(x = x_prop_lps, beta = beta_schedule[K]), target_args))
        lsaq_c2p <- lpsampler_q(x = x_prop_lps)
        lsaq_p2c <- lpsampler_q(x = x_curr_lps)
        lsa_step <- mh_step(x_curr = x_curr_lps, x_prop = x_prop_lps,
                            l_curr = l_curr_lps, l_prop = l_prop_lps,
                            lq_c2p = lsaq_c2p, lq_p2c = lsaq_p2c,
                            do_checks = FALSE)

        mh_acc[c, s, K] <- lsa_step$accepted
        ww_n <- i + s
        ww_p <- ww_n - 1
        x[ww_n, k_i, ] <- lsa_step$x_next
        l_x[ww_n, k_i] <- lsa_step$l_next
        y[ww_p, k_i, ] <- lsa_step$x_prop
        l_y[ww_p, k_i] <- lsa_step$l_prop
        delta_l[ww_p, k_i] <- lsa_step$delta_l
      }

    }else{

      mh_level_args <- list(x_0 = x[i, k_i , ],
                            l_0 = l_x[i, k_i],
                            beta = beta_schedule[K],
                            mh_sampler = sampler_list[[K]],
                            other_sampler_args = sampler_args_list[[K]],
                            lq_mh = lq_list[[K]],
                            other_lq_args = lq_args_list[[K]],
                            S = Within_Moves, burn = 0, silent = TRUE, d = d)
      mh_moves <- do.call(mh_sampler_chain, c(mh_level_args, l_target, target_args))

      mh_acc[c, , K] <- mh_moves$acc
      x[within_window_next, k_i, ] <- mh_moves$x
      l_x[within_window_next, k_i] <- mh_moves$l_x
      y[within_window_prop, k_i, ] <- mh_moves$y
      l_y[within_window_prop, k_i] <- mh_moves$l_y
      delta_l[within_window_prop, k_i] <- mh_moves$delta_l

    }


    # Copy parameters for next cycle
    k_indexes[c+1, 1, ] <- k_indexes[c, Temp_Moves + 1, ]
    beta_indexes[c+1, 1, ] <- beta_indexes[c, Temp_Moves + 1, ]

  }
  global_times[2] <- Sys.time()

  #--- Post processing and result -------------

  swap_acc_rates <- apply(swap_acc, 3, function(swaps) mean(swaps[swaps>=0]))
  mh_acc_rates <- colMeans(mh_acc, dims = 2)

  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Swap Acceptance Rates:", round(swap_acc_rates,3)), sep = "\n")
    cat(paste("MH Acceptance Rates:", round(mh_acc_rates,3)), sep = "\n")
  }

  if(burn_cycles == -1){
    global_times[3] <- Sys.time()
    return(list(x = x, l_x = l_x, y = y, l_y = l_y, delta_l = delta_l, u_lps = u_lps,
                k_indexes = k_indexes[-(Cycles+1), , ],
                beta_indexes = beta_indexes[-(Cycles+1), , ],
                swap_acc = swap_acc, swap_acc_rates = swap_acc_rates,
                mh_acc = mh_acc, mh_acc_rates = mh_acc_rates,
                cycle_times = cycle_times, global_times = global_times))
  }

  burn_window <- seq(1,burn_cycles*cycle_length + 1)
  x_r <- x[-burn_window, , , drop = FALSE]
  l_x_r <- l_x[-burn_window, ]
  y_r <- y[-burn_window, , , drop = FALSE]
  l_y_r <- l_y[-burn_window, ]
  delta_l_r <- delta_l[-burn_window, ]
  u_lps_r <- u_lps[-burn_window, ]

  if(burn_cycles == 0){
    mh_acc_r <- mh_acc
    swap_acc_r <- swap_acc
    b_r <- beta_indexes[-c(1,Cycles+1), , ]
    k_r <- k_indexes[-c(1,Cycles+1), , ]
  }else{
    burn_window <- 1:burn_cycles
    mh_acc_r <- mh_acc[-burn_window, , ]
    swap_acc_r <- swap_acc[-burn_window, ]
    b_r <- beta_indexes[-c(burn_window,Cycles + 1), , ]
    k_r <- k_indexes[-c(burn_window,Cycles + 1), , ]
  }

  global_times[3] <- Sys.time()
  return(list(x = x_r, l_x = l_x_r, y = y_r, l_y = l_y_r, delta_l = delta_l_r, u_lps = u_lps_r,
              k_indexes = k_r, beta_indexes = b_r,
              swap_acc = swap_acc_r, swap_acc_rates = swap_acc_rates,
              mh_acc = mh_acc_r, mh_acc_rates = mh_acc_rates,
              cycle_times = cycle_times, global_times = global_times))

}

#' @export
ALPS_rwm_leaner_chain_list <- function(ltemp_target, ..., HAT_info,
                                       beta_schedule, swap_type = "deo", quanta_levels = NULL,
                                       scale = 1, Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                                       x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, jump_p = 0.7,
                                       k_0 = NULL, custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                                       quanta_mode_info = NULL,
                                       quanta_pass_prev_mod_assign = TRUE,
                                       silent = FALSE){

  #--- HAT use -------------
  l_target <- lHAT_target_cpp
  target_args <- c(list(ltemp_target = ltemp_target),
                   "HAT_info" = rlang::dots_list(HAT_info),
                   rlang::dots_list(...))


  #--- Preparation -------------

  # General PT parameters
  start_time <- Sys.time()
  global_times <- rep(start_time, 3)
  K <- length(beta_schedule)
  beta_max <- beta_schedule[K]
  odd_indices <- seq(1, K, by = 2)
  even_indices <- seq(2, K, by = 2)


  # Within Level Moves Dispatcher
  make_within_level_moves <- function(x_w, l_w, beta_w, scale_w, sampler_w, u_w){
    if(u_w){
      # Leap-Point Sampler
      lps_args_list <- list(x_0 = x_w, l_0 = l_w, beta = beta_w,
                            mh_sampler = lpsampler_cpp,
                            other_sampler_args = c("beta_max" = beta_w,
                                                   HAT_info[c("w","modes","L")]),
                            lq_mh = lps_q_cpp,
                            other_lq_args = c("beta_max" = beta_w,
                                              HAT_info[c("w","modes","L_inv","ldet_L_inv")]),
                            S = Within_Moves, d = d, burn = 0, silent = TRUE) |>
        c("l_target" = l_target, target_args)
      aux <- do.call(mh_sampler_leaner_chain_list, lps_args_list)
    }else{
      # Regular RWM
      rwm_args_list <- list(x_0 = x_w, l_0 = l_w, beta = beta_w,
                            scale = scale_w, custom_rw_sampler = sampler_w,
                            S = Within_Moves, d = d, burn = 0, silent = TRUE) |>
        c("l_target" = l_target, target_args)
      aux <- do.call(rwm_sampler_leaner_chain_list, rwm_args_list)
    }
    return(aux)
  }


  # QuanTA levels checking
  if(is.null(quanta_levels)){
    quanta_levels <- K-1
  }else{
    if(!all(quanta_levels %in% seq(1,K-1))){
      stop(paste("Invalid quanta_levels.",
                 "These must be the base indexes for which QuanTA swaps will be attempted,",
                 "hence a vector of integers between 1 and K-1,",
                 "where K is the number of temperatures used.",
                 sep = "\n"))
    }
    quanta_mode_info <- quanta_mode_info %||% HAT_info
  }

  # Dimension and Scales
  if(is.list(scale)){
    if(length(scale) != K){
      stop("As a list, length of scale must must equal that of beta_schedule")
    }
    scale_list <- scale
  }else{
    stopifnot(is.numeric(scale))
    scale_list <- lapply(beta_schedule, function(beta) scale/beta)
  }
  d <- d %||% ifelse(is.matrix(scale_list[[1]]), nrow(scale_list[[1]]), length(scale_list[[1]]))
  if(d > 1 && !is.matrix(scale_list[[1]])){
    scale_list <- lapply(scale_list, function(scale) diag(scale, d))
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
  }

  # Checking for valid sample sizes
  stopifnot(Cycles >= 1)
  stopifnot(-1 <= burn_cycles && burn_cycles < Cycles)
  stopifnot(Temp_Moves >= 1)
  stopifnot(Within_Moves >= 1)

  # If the user didn't we define proposal sampler(s) as indep. normals
  if(is.list(custom_rw_sampler)){
    if(length(custom_rw_sampler) != K){
      stop("As a list, custom_rw_sampler must have the same length of beta_schedule")
    }
    sampler_list <- custom_rw_sampler
  }else{
    sampler <- custom_rw_sampler %||%
      ifelse(d == 1,
             function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
             function(x, scale){ rmvtnorm_temp_cpp(n = 1, mu = x, sigma = scale) })
    sampler_list <- rep(list(sampler),K)
  }

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- matrix(data = runif(d*K, min = -x_0_u, max = x_0_u), nrow = d, ncol = K)
  }else{
    stopifnot(is.numeric(x_0))
    stopifnot(nrow(x_0) == d && ncol(x_0) == K)
  }
  k_0 <- k_0 %||% 1:K
  b_0 <- beta_schedule[k_0]
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(m in 1:K){
      l_0[m] <- do.call(l_target, c(list(x = x_0[ , m], beta = b_0[m]), target_args))
    }
  }


  # Preallocate containers
  cycle_length <- Within_Moves + Temp_Moves
  x <- lapply(1:K, function(m) matrix(NA_real_, nrow = d, ncol = Cycles*cycle_length + 1))
  k_indexes <- matrix(NA_integer_, nrow = Cycles*Temp_Moves + 1, ncol = K)
  beta_indexes <- matrix(NA_real_, nrow = Cycles*Temp_Moves + 1, ncol = K)
  swap_acc <- array(-1, dim = c(K-1, Temp_Moves, Cycles))
  rwm_acc <- lapply(1:K, function(k) matrix(NA, nrow = Within_Moves, ncol = Cycles))

  #--- Algorithm -------------

  # Initialize
  for(m in 1:K){
    x[[m]][ , 1] <- x_0[ , m]
  }
  k_indexes[1, ] <- k_0
  beta_indexes[1, ] <- b_0
  l_x <- l_0


  # Run Cycles
  i_cycle <- 1
  j_cycle <- 1
  sequential_plan <- is(future::plan(), "sequential") # possible parallelism
  for(c in 1:Cycles){

    if(!silent & isTRUE(c %% floor(Cycles*0.05) == 0)){
      cat(paste0("Avance: ", round(100*c/Cycles), "%"), sep = "\n")
      cat(format(Sys.time(), usetz=TRUE), sep = "\n")
    }

    # Within Level Exploration (possible parallelism)
    if(sequential_plan){
      within_level_moves <- mapply(
        FUN = make_within_level_moves,
        x_w = lapply(x, function(x_m) x_m[, i_cycle]),
        l_w = l_x,
        beta_w = beta_indexes[j_cycle, ],
        scale_w = scale_list[k_indexes[j_cycle, ]],
        sampler_w = sampler_list[k_indexes[j_cycle, ]],
        u_w = (runif(1) <= jump_p) & (k_indexes[j_cycle, ] == K),
        SIMPLIFY = FALSE)
    }else{
      within_level_moves <- future.apply::future_mapply(
        FUN = make_within_level_moves,
        x_w = lapply(x, function(x_m) x_m[, i_cycle]),
        l_w = l_x,
        beta_w = beta_indexes[j_cycle, ],
        scale_w = scale_list[k_indexes[j_cycle, ]],
        sampler_w = sampler_list[k_indexes[j_cycle, ]],
        u_w = u_lps[c] & (k_indexes[j_cycle, ] == K),
        SIMPLIFY = FALSE,
        future.seed = TRUE,
        future.globals = c("Within_Moves","d","HAT_info"))
    }
    for(m in 1:K){
      x[[m]][ , i_cycle + 1:Within_Moves] <- t(within_level_moves[[m]]$x)
      l_x[m] <- within_level_moves[[m]]$l_x_curr
      rwm_acc[[ k_indexes[j_cycle, m] ]][ , c] <- within_level_moves[[m]]$acc
    }


    i_cycle <- i_cycle + Within_Moves

    # Temperature Swaps
    for(t in 1:Temp_Moves){

      swap_arguments <- list(x_curr = lapply(x, function(x_m) x_m[, i_cycle]),
                             l_curr = l_x,
                             beta_curr = beta_indexes[j_cycle, ],
                             k_curr = k_indexes[j_cycle, ],
                             type = swap_type, j_deo = ifelse(Temp_Moves == 1, c, t),
                             l_target, target_args, quanta_levels = quanta_levels,
                             mode_info = quanta_mode_info, K = K, d = d,
                             odd_indices = odd_indices, even_indices = even_indices,
                             pass_mod_assignment = quanta_pass_prev_mod_assign)
      swap_moves <- do.call(alps_swap_move_list, swap_arguments)

      swap_acc[, t, c] <- swap_moves$acc
      l_x <- swap_moves$l_next

      i_cycle <- i_cycle + 1
      for(m in 1:K){
        x[[m]][, i_cycle] <- swap_moves$x_next[[m]]
      }

      j_cycle <- j_cycle + 1
      k_indexes[j_cycle, ] <- swap_moves$k_next
      beta_indexes[j_cycle, ] <- swap_moves$beta_next

    }

  }
  global_times[2] <- Sys.time()

  #--- Post processing and result -------------

  swap_acc_rates <- apply(swap_acc, 1, function(swaps) mean(swaps[swaps>=0]) )
  rwm_acc_rates <- vapply(rwm_acc, mean, numeric(1))

  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Swap Acceptance Rates:", round(swap_acc_rates,3)), sep = "\n")
    cat(paste("RWM Acceptance Rates:", round(rwm_acc_rates,3)), sep = "\n")
  }

  if(burn_cycles == -1){
    global_times[3] <- Sys.time()
    return(mget(c("x","k_indexes","beta_indexes",
                  "swap_acc","swap_acc_rates",
                  "rwm_acc","rwm_acc_rates",
                  "global_times")))
  }

  burn_window <- seq(1, burn_cycles + 1)
  burn_window_k <- seq(1, burn_cycles*Temp_Moves + 1)
  global_times[3] <- Sys.time()
  return(list(x = lapply(x,function(x_m) x_m[,-burn_window]),
              k_indexes = k_indexes[-burn_window_k],
              beta_indexes =  beta_indexes[-burn_window_k],
              swap_acc = swap_acc[,,-burn_window],
              swap_acc_rates = swap_acc_rates,
              rwm_acc = lapply(rwm_acc, function(moves) moves[,-burn_window]),
              rwm_acc_rates = rwm_acc_rates,
              global_times = global_times))

}

