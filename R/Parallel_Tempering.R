
PT_rwm_chain <- function(l_target, ..., beta_schedule, swap_type = "deo",
                         scale = 1, Temp_Moves = 1000, Within_Moves = 10, burn_cycles = 0,
                         x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                         custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                         silent = FALSE){

#--- Preparation -------------

  # General PT parameters
  K <- length(beta_schedule)
  stopifnot(swap_type %in% c("deo","seo","naive"))
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
  stopifnot(Temp_Moves >= 1)
  stopifnot(Within_Moves >= 1)
  stopifnot(0 <= burn_cycles && burn_cycles < Temp_Moves)

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
  cycle_length <- Within_Moves + 1
  S_Tot <- Temp_Moves*(cycle_length)
  k_indexes <- matrix(nrow = Temp_Moves + 1, ncol = K)
  beta_indexes <- matrix(nrow = Temp_Moves + 1, ncol = K)
  x <- array(dim = c(S_Tot + 1, K, d), dimnames = list(NULL, NULL, target_names))
  l <- matrix(nrow = S_Tot + 1, ncol = K)
  swap_acc <- matrix(nrow = Temp_Moves, ncol = K)
  rwm_acc <- array(dim = c(Temp_Moves, Within_Moves, K))

#--- Algorithm -------------

  # Initialize
  k_indexes[1, ] <- 1:K
  beta_indexes[1, ] <- beta_schedule
  x[1, , ] <- x_0
  l[1, ] <- l_0

  # Run iterations
  for(c in 1:Temp_Moves){

    if(c %% floor(Temp_Moves*0.05) == 0){
      cat(paste0("Avance: ",round(100*c/Temp_Moves),"%"),sep = "\n")
      print(Sys.time())
    }
    # Cycle index
    i <- (c-1)*(cycle_length) + 2

    # Temperature Swap
    swap_move <- temp_swap_move(type = swap_type, c = c,
                                beta_curr = beta_indexes[c, ],
                                k_curr = k_indexes[c, ],
                                x_curr = x[i-1, , , drop = FALSE],
                                l_curr = l[i-1, ],
                                l_target, target_args,
                                K = K, d = d)
    swap_acc[c, ] <- swap_move$acc
    k_indexes[c+1, ] <- swap_move$k_next
    beta_indexes[c+1, ] <- swap_move$beta_next
    x[i, , ] <- swap_move$x_next
    l[i, ] <- swap_move$l_next

    # Within temperature moves
    for(k in 1:K){
      temperature_arguments <- c(list(x_0 = x[i, k , ],
                                      l_0 = l[i, k],
                                      beta = beta_indexes[c+1, k],
                                      custom_rw_sampler = sampler_list[[ swap_move$k_next[k] ]],
                                      scale = scale_list[[ swap_move$k_next[k] ]],
                                      S = Within_Moves, burn = 0, silent = TRUE),
                                 l_target = l_target, target_args)
      rwm_moves <- do.call(rwm_sampler_chain, temperature_arguments)
      rwm_acc[c, , swap_move$k_next[k]] <- rwm_moves$acc
      x[i + 1:Within_Moves, k, ] <- rwm_moves$x
      l[i + 1:Within_Moves, k] <- rwm_moves$l
    }

  }

#--- Post processing and result -------------

  swap_acc_rates <- colMeans(swap_acc, na.rm = TRUE)
  rwm_acc_rates <- colMeans(rwm_acc, na.rm = TRUE, dims = 2)

  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Swap Acceptance Rates:", round(swap_acc_rates,3)), sep = "\n")
    cat(paste("RWM Acceptance Rates:", round(rwm_acc_rates,3)), sep = "\n")
  }

  x_r <- x[-seq(1,burn_cycles*cycle_length + 1), , , drop = FALSE]
  l_r <- l[-seq(1,burn_cycles*cycle_length + 1), ]

  if(burn_cycles == 0){
    rwm_acc_r <- rwm_acc
    swap_acc_r <- swap_acc
    b_r <- beta_indexes[-1, ]
    k_r <- k_indexes[-1, ]
  }else{
    rwm_acc_r <- rwm_acc[-seq(1,burn_cycles), , ]
    swap_acc_r <- swap_acc[-seq(1,burn_cycles), ]
    b_r <- beta_indexes[-seq(1,1+burn_cycles), ]
    k_r <- k_indexes[-seq(1,1+burn_cycles), ]
  }

  return(list(x = x_r, l = l_r, k_indexes = k_r, beta_indexes = b_r,
              swap_acc = swap_acc_r, swap_acc_rates = swap_acc_rates,
              rwm_acc = rwm_acc_r, rwm_acc_rates = rwm_acc_rates))

}
