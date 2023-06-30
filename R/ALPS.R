
#' @export
ALPS_rwm_chain <- function(ltemp_target, ..., HAT = TRUE, HAT_info,
                           beta_schedule, swap_type = "deo", quanta_levels = NULL,
                           scale = 1, Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                           x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, jump_p = 0.7,
                           custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                           quanta_mode_info = NULL, silent = FALSE){

  #--- HAT use -------------
  if(HAT){
    l_target <- lHAT_target
    target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                     rlang::dots_list(...))
  }else{
    l_target <- ltemp_target
    target_args <- rlang::dots_list(...)
  }

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
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(k in seq_along(beta_schedule)){
      temporal <- do.call(l_target, c(list(x = x_0[k, ], beta = beta_schedule[k]), target_args))
      l_0[k] <- temporal
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
  k_indexes[1, 1, ] <- 1:K
  beta_indexes[1, 1, ] <- beta_schedule
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
                                  custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                                  quanta_mode_info = NULL, silent = FALSE){

  #--- HAT use -------------
  if(HAT){
    l_target <- lHAT_target
    target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                     rlang::dots_list(...))
  }else{
    l_target <- ltemp_target
    target_args <- rlang::dots_list(...)
  }

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
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(k in seq_along(beta_schedule)){
      temporal <- do.call(l_target, c(list(x = x_0[k, ], beta = beta_schedule[k]), target_args))
      l_0[k] <- temporal
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
  k_indexes[1, 1, ] <- 1:K
  beta_indexes[1, 1, ] <- beta_schedule
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
      rwm_moves <- do.call(rwm_sampler_chain, c(rwm_level_args, l_target, target_args))

      rwm_acc[c, , k] <- rwm_moves$acc
      x[within_window_next, k_i, ] <- rwm_moves$x
      l_x[k_i] <- rwm_moves$l_x[Within_Moves]

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
                          quanta_mode_info = NULL, silent = FALSE, target_names = NULL){

#--- HAT use -------------
  if(HAT){
    l_target <- lHAT_target
    target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                     rlang::dots_list(...))
  }else{
    l_target <- ltemp_target
    target_args <- rlang::dots_list(...)
  }

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
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(k in seq_along(beta_schedule)){
      temporal <- do.call(l_target, c(list(x = x_0[k, ], beta = beta_schedule[k]), target_args))
      l_0[k] <- temporal
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
  k_indexes[1, 1, ] <- 1:K
  beta_indexes[1, 1, ] <- beta_schedule
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
ALPS_rwm_leaner_chain_list <- function(ltemp_target, ..., HAT = TRUE, HAT_info,
                                       beta_schedule, swap_type = "deo", quanta_levels = NULL,
                                       scale = 1, Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                                       x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, jump_p = 0.7,
                                       k_0 = NULL, custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                                       quanta_mode_info = NULL, silent = FALSE){

  #--- HAT use -------------
  if(HAT){
    l_target <- lHAT_target
    target_args <- c(list(HAT_info = HAT_info, ltemp_target = ltemp_target),
                     rlang::dots_list(...))
  }else{
    l_target <- ltemp_target
    target_args <- rlang::dots_list(...)
  }

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

  # Within Level Moves Dispatcher
  make_within_level_moves <- function(x_w, l_w, beta_w, scale_w, sampler_w, u_w,
                                      l_target_w, other_target_args, mh_sampler_w, lq_mh_w,
                                      S_w, d_w){
    if(u_w){
      lps_args_list <- list(x_0 = x_w, l_0 = l_w, l_target = l_target_w,
                            mh_sampler = mh_sampler_w, lq_mh = lq_mh_w,
                            S = S_w, d = d_w, burn = 0, silence = TRUE) |>
        c(other_target_args)
      return(do.call(mh_sampler_leaner_chain,lps_args_list))
    }
    rwm_args_list <- list(x_0 = x_w, l_0 = l_w, l_target = l_target_w, beta = beta_w,
                          scale = scale_w, custom_rw_sampler = sampler_w,
                          S = S_w, d = d_w, burn = 0, silent = TRUE) |>
      c(other_target_args)
    return(do.call(rwm_sampler_chain, rwm_args_list))

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
    x_0 <- matrix(data = runif(d*K, min = -x_0_u, max = x_0_u), nrow = d, ncol = K)
  }else{
    stopifnot(is.numeric(x_0))
    stopifnot(nrow(x_0) == d && ncol(x_0) == K)
  }
  if(is.null(k_0)){
    k_0 <- 1:K
    b_0 <- beta_schedule
  }else{
    b_0 <- beta_schedule[k_0]
  }
  if(is.null(l_0)){
    # More verbose but attempts to use mapply have failed
    # (asplit for the matrix rows may also fail because it keeps them as arrays)
    l_0 <- numeric(K)
    for(m in 1:K){
      l_0[m] <- do.call(l_target, c(list(x = x_0[ , k_0[m]], beta = b_0[m]), target_args))
    }
  }


  # Preallocate containers
  cycle_length <- Temp_Moves + Within_Moves
  x <- lapply(1:(Cycles + 1), function(c){
    lapply(1:K, function(k) matrix(NA_real_, cycle_length, d))
  })
  k_indexes <- lapply(1:(Cycles + 1), function(c) matrix(NA_integer_, nrow = Temp_Moves, ncol = K))
  beta_indexes <- lapply(1:(Cycles + 1), function(c) matrix(NA_real_, nrow = Temp_Moves, ncol = K))
  swap_acc <- array(-1, dim = c(K-1, Temp_Moves, Cycles))
  rwm_acc <- lapply(1:K, function(k) matrix(NA, nrow = Within_Moves, ncol = Cycles))

  #--- Algorithm -------------

  # Initialize
  for(m in 1:K){
    x[[1]][[m]][cycle_length, ] <- x_0[ , m]
  }
  k_indexes[[1]][Temp_Moves, ] <- k_0
  beta_indexes[[1]][Temp_Moves, ] <- b_0
  l_x <- l_0

  # Determine LPS cycles at coldest level
  u_lps <- runif(Cycles) <= jump_p


  # Run Cycles
  for(c in 1:Cycles){

    i_c <- c+1
    if(!silent & isTRUE(c %% floor(Cycles*0.05) == 0)){
      cat(paste0("Avance: ", round(100*c/Cycles), "%"), sep = "\n")
      cat(format(Sys.time(), usetz=TRUE), sep = "\n")
    }

    # Within Level Exploration
    # Iteration containers (for possible parallelism)
    x_c <- lapply(x[[c]],function(x_m) x_m[cycle_length, ])
    k_c <- k_indexes[[c]][Temp_Moves, ]
    u_c <- u_lps[c] & k_c == K
    beta_c <- beta_indexes[[c]][Temp_Moves, ]
    scale_c <- scale_list[k_c]
    sampler_c <- sampler_list[k_c]
    within_level_moves <- mapply(FUN = make_within_level_moves,
                                 x_w = x_c, l_w = l_x, beta_w = beta_c,
                                 scale_w = scale_c,  sampler_w = sampler_c,
                                 u_w = u_c,
                                 MoreArgs = list(l_target_w = l_target,
                                                 other_target_args = target_args,
                                                 mh_sampler_w = lpsampler, lq_mh_w = lpsampler_q,
                                                 S_w = Within_Moves, d_w = d),
                                 SIMPLIFY = FALSE)
    for(m in 1:K){
      x[[i_c]][[m]][1:Within_Moves,] <- within_level_moves[[m]]$x
      l_x[m] <- within_level_moves[[m]]$l_x_curr
      rwm_acc[[k_c[m]]][,c] <- within_level_moves[[m]]$acc
    }

    # Temperature Swaps
    swap_arguments <- list(type = swap_type, j_deo = ifelse(Temp_Moves == 1, c, 1),
                           quanta_levels = quanta_levels, mode_info = quanta_mode_info,
                           beta_curr = beta_c, k_curr = k_c,
                           x_curr = x_c,
                           l_curr = l_x,
                           l_target, target_args,
                           K = K,
                           odd_indices = odd_indices,
                           even_indices = even_indices,
                           d = d)
    swap_move <- do.call(alps_swap_move_list, swap_arguments)

    swap_acc[, 1, c] <- swap_move$acc
    k_indexes[[i_c]][1, ] <- swap_move$k_next
    beta_indexes[[i_c]][1, ] <- swap_move$beta_next
    for(m in 1:K){
      x[[i_c]][[m]][Within_Moves + 1, ] <- swap_move$x_next[[m]]
    }
    l_x <- swap_move$l_next
    if(Temp_Moves > 1){
      for(j in 2:Temp_Moves){

        swap_arguments <- list(type = swap_type, j_deo = j,
                               quanta_levels = quanta_levels, mode_info = quanta_mode_info,
                               beta_curr = beta_indexes[[i_c]][j-1, ],
                               k_curr = k_indexes[[i_c]][j-1, ],
                               x_curr = lapply(x[[i_c]], function(x_m) x_m[Within_Moves + j - 1, ]),
                               l_curr = l_x,
                               l_target, target_args,
                               K = K,
                               odd_indices = odd_indices,
                               even_indices = even_indices,
                               d = d)
        swap_move <- do.call(alps_swap_move_list, swap_arguments)

        swap_acc[, j, c] <- swap_move$acc
        k_indexes[[i_c]][j, ] <- swap_move$k_next
        beta_indexes[[i_c]][j, ] <- swap_move$beta_next
        for(m in 1:K){
          x[[i_c]][[m]][Within_Moves + j, ] <- swap_move$x_next[[m]]
        }
        l_x <- swap_move$l_next

      }
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
  global_times[3] <- Sys.time()
  return(list(x = x[-burn_window],
              k_indexes = k_indexes[-burn_window],
              beta_indexes =  beta_indexes[-burn_window],
              swap_acc = swap_acc[,,-burn_window],
              swap_acc_rates = swap_acc_rates,
              rwm_acc = lapply(rwm_acc, function(moves) moves[,-burn_window]),
              rwm_acc_rates = rwm_acc_rates,
              global_times = global_times))

