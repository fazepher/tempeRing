
ST_rwm_chain <- function(l_target, ..., beta_schedule, g_schedule = NULL,
                         scale = 1, Temp_Moves = 1000, Within_Moves = 10, burn_cycles = 0,
                         x_0 = NULL, x_0_u = 2, k_0 = NULL, l_0 = NULL, seed = NULL,
                         custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                         silent = FALSE){

#--- Preparation -------------

  # General ST parameters
  K <- length(beta_schedule)
  g_schedule <- g_schedule %||% rep(0, K)
  if(length(g_schedule) != K){
    stop("Length of g_schedule must match that of beta_schedule")
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
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }
  if(is.null(k_0) && !is.null(l_0)){
    stop("When l_0 is provided, k_0 must also be present")
  }


  # Preallocate containers
  cycle_length <- Within_Moves + 1
  S_Tot <- Temp_Moves*(cycle_length)
  k <- numeric(Temp_Moves + 1)
  x <- matrix(nrow = S_Tot + 1, ncol = d, dimnames = list(NULL, target_names))
  l <- numeric(S_Tot + 1)
  swap_acc <- logical(Temp_Moves)
  rwm_acc <- matrix(nrow = Temp_Moves, ncol = Within_Moves)

#--- Algorithm -------------

  # Initialize
  k[1] <- k_0 %||% sample(1:K, 1)
  x[1, ] <- x_0
  l[1] <- l_0 %||% l_target(x_0, beta_schedule[k[1]], ...)

  # Run iterations
  for(c in 1:Temp_Moves){

    # Cycle index
    i <- (c-1)*(cycle_length) + 2

    # Temperature Swap
    temp_move <- st_temp_step(k_curr = k[c],
                              x_curr = x[i-1, ],
                              l_curr = l[i-1],
                              l_target = l_target, ...,
                              beta_schedule = beta_schedule, g_schedule = g_schedule, K = K)
    swap_acc[c] <- temp_move$acc
    k[c + 1] <- temp_move$k_next
    x[i, ] <- x[i-1, ] # We always remain at the same X state in ST
    l[i] <- temp_move$l_next

    # Within temperature moves
    rwm_moves <- rwm_sampler_chain(x_0 = x[i, ],
                                   l_0 = l[i],
                                   beta = beta_schedule[ temp_move$k_next ],
                                   custom_rw_sampler = sampler_list[[ temp_move$k_next ]],
                                   scale = scale_list[[ temp_move$k_next ]],
                                   l_target = l_target, ...,
                                   S = Within_Moves, burn = 0, silent = TRUE)
    rwm_acc[c, ] <- rwm_moves$acc
    x[i + 1:Within_Moves, ] <- rwm_moves$x
    l[i + 1:Within_Moves] <- rwm_moves$l

  }

#--- Post processing and result -------------

  swap_acc_rates <- vapply(1:K, function(level) swap_acc[which(k[1:Temp_Moves] == level)] |> mean(),
                           numeric(1))
  rwm_acc_rates <- vapply(1:K, function(level) rwm_acc[which(k[-1] == level), ] |> mean(),
                          numeric(1))

  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Swap Acceptance Rates:", round(swap_acc_rates,3)), sep = "\n")
    cat(paste("RWM Acceptance Rates:", round(rwm_acc_rates,3)), sep = "\n")
  }

  x_r <- x[-seq(1,burn_cycles*cycle_length + 1), ]
  l_r <- l[-seq(1,burn_cycles*cycle_length + 1)]

  if(burn_cycles == 0){
    rwm_acc_r <- rwm_acc
    swap_acc_r <- swap_acc
    k_r <- k[-Temp_Moves]
  }else{
    rwm_acc_r <- rwm_acc[-seq(1,burn_cycles), ]
    swap_acc_r <- swap_acc[-seq(1,burn_cycles)]
    k_r <- k[-seq(1,1+burn_cycles)]
  }

  return(list(x = x_r, k = k_r, l = l_r,
              swap_acc = swap_acc_r, swap_acc_rates = swap_acc_rates,
              rwm_acc = rwm_acc_r, rwm_acc_rates = rwm_acc_rates))

}
