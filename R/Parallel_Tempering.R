
#' @export
PT_rwm_chain <- function(l_target, ..., beta_schedule, swap_type = "deo",
                         scale = 1, Cycles = 1000, Temp_Moves = 5, Within_Moves = 5, burn_cycles = 0,
                         x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                         custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                         quanta = FALSE, mode_info = NULL,
                         silent = FALSE){

#--- Preparation -------------

  # General PT parameters
  start_time <- Sys.time()
  global_times <- rep(start_time, 3)
  K <- length(beta_schedule)
  odd_indices <- seq(1, K, by = 2)
  even_indices <- seq(2, K, by = 2)

  # Dimension and Scales
  if(is.list(scale)){
    if(length(scale) != K){
      stop("As a list, length of scale must must match that of beta_schedule")
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
      stop("As a list, custom_rw_sampler must have the same length as beta_schedule")
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
  target_args <- rlang::dots_list(...)
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
  l <- matrix(nrow = S_Tot + 1, ncol = K)
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
  l[1, ] <- l_0

  # Run iterations
  for(c in 1:Cycles){

    cycle_times[c] <- Sys.time()
    if(isTRUE(c %% floor(Cycles*0.05) == 0)){
      cat(paste0("Avance: ",round(100*c/Cycles),"%"),sep = "\n")
      print(Sys.time())
    }

    # Current cycle position index
    i <- (c-1)*(cycle_length) + 2

    # Temperature Swaps
    for(j in 1:Temp_Moves){

      swap_arguments <- list(type = swap_type, j_deo = ifelse(Temp_Moves == 1, c, j),
                             quanta = quanta, mode_info = mode_info,
                             beta_curr = beta_indexes[c, j, ],
                             k_curr = k_indexes[c, j,  ],
                             x_curr = x[(i-2)+j, , , drop = FALSE],
                             l_curr = l[(i-2)+j, ],
                             l_target, target_args,
                             K = K,
                             odd_indices = odd_indices,
                             even_indices = even_indices,
                             d = d)
      swap_move <- do.call(temp_swap_move, swap_arguments)

      swap_acc[c, j, ] <- swap_move$acc
      k_indexes[c, j+1, ] <- swap_move$k_next
      beta_indexes[c, j+1, ] <- swap_move$beta_next

      x[(i-1) + j, , ] <- swap_move$x_next
      l[(i-1) + j, ] <- swap_move$l_next

    }

    # Within temperature moves

    # Update current cycle position index
    i <- i + Temp_Moves - 1
    for(k in 1:K){

      rwm_level_args <- list(x_0 = x[i, k , ],
                             l_0 = l[i, k],
                             beta = beta_indexes[c, Temp_Moves + 1, k],
                             custom_rw_sampler = sampler_list[[ k_indexes[c, Temp_Moves + 1, k] ]],
                             scale = scale_list[[ k_indexes[c, Temp_Moves + 1, k] ]],
                             S = Within_Moves, burn = 0, silent = TRUE)
      rwm_moves <- do.call(rwm_sampler_chain, c(rwm_level_args, l_target, target_args))

      rwm_acc[c, , k_indexes[c, Temp_Moves + 1, k]] <- rwm_moves$acc
      x[i + 1:Within_Moves, k, ] <- rwm_moves$x
      l[i + 1:Within_Moves, k] <- rwm_moves$l

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
    return(list(x = x, l = l,
                k_indexes = k_indexes[-(Cycles+1), , ],
                beta_indexes = beta_indexes[-(Cycles+1), , ],
                swap_acc = swap_acc, swap_acc_rates = swap_acc_rates,
                rwm_acc = rwm_acc, rwm_acc_rates = rwm_acc_rates,
                cycle_times = cycle_times, global_times = global_times))
  }

  x_r <- x[-seq(1,burn_cycles*cycle_length + 1), , , drop = FALSE]
  l_r <- l[-seq(1,burn_cycles*cycle_length + 1), ]

  if(burn_cycles == 0){
    rwm_acc_r <- rwm_acc
    swap_acc_r <- swap_acc
    b_r <- beta_indexes[-c(1,Cycles+1), , ]
    k_r <- k_indexes[-c(1,Cycles+1), , ]
  }else{
    rwm_acc_r <- rwm_acc[-seq(1,burn_cycles), , ]
    swap_acc_r <- swap_acc[-seq(1,burn_cycles), ]
    b_r <- beta_indexes[-c(seq(1,burn_cycles),Cycles + 1), , ]
    k_r <- k_indexes[-c(seq(1,burn_cycles),Cycles + 1), , ]
  }

  global_times[3] <- Sys.time()
  return(list(x = x_r, l = l_r, k_indexes = k_r, beta_indexes = b_r,
              swap_acc = swap_acc_r, swap_acc_rates = swap_acc_rates,
              rwm_acc = rwm_acc_r, rwm_acc_rates = rwm_acc_rates,
              cycle_times = cycle_times, global_times = global_times))

}
