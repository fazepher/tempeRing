
rwm_sampler_chain <- function(l_target, ..., scale = 1, S = 1000, burn = 0,
                              x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                              custom_rw_sampler = NULL, more_sampler_args = NULL,
                              d = NULL, silent = FALSE){

#--- Preparation -------------

  # Define the dimension and sampler scale
  stopifnot(is.numeric(scale))
  d <- d %||% ifelse(is.matrix(scale), nrow(scale), length(scale))
  if(d > 1 && !is.matrix(scale)){
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
    scale <- diag(scale)
  }

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(-1 <= burn && burn < S)

  # If the user didn't, we define the proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ rmvtnorm(n = 1, mu = x, sigma = scale) })

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

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  y <- matrix(nrow = S, ncol = d)
  l_x <- numeric(S + 1)
  l_y <- numeric(S)
  delta_l <- numeric(S)
  acc <- logical(S)

#--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x[1] <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step(x[s, ], l_x[s], l_target, ...,
                                         sampler = sampler,
                                         sampler_args = c(list(scale), more_sampler_args),
                                         do_checks = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_x[s+1] <- rwm_step$l_next
    y[s, ] <- rwm_step$x_prop
    l_y[s] <- rwm_step$l_prop
    delta_l[s] <- rwm_step$delta_l
    acc[s] <- rwm_step$accepted
  }

#--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(mget(c("x", "l_x", "acc", "acc_rate", "y", "l_y", "delta_l")))
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x = l_x[-1],
                acc = acc, acc_rate = acc_rate,
                y = y, l_y = l_y, delta_l = delta_l))
  }

  bwx <- 1:(burn+1)
  bwy <- 1:burn
  return(list(x = x[-bwx, ], l_x = l_x[-bwx],
              acc = acc[-bwy], acc_rate = acc_rate,
              y = y[-bwy, ], l_y = l_y[-bwy], delta_l = delta_l[-bwy]))

}

mh_sampler_chain <- function(l_target, ...,
                             mh_sampler, other_sampler_args = NULL,
                             lq_mh, other_lq_args = NULL,
                             S = 1000, burn = 0, d = 1,
                             x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL, silent = FALSE){

  #--- Preparation -------------

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(0 <= burn && burn < S)

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

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d)
  y <- matrix(nrow = S, ncol = d)
  l_x <- numeric(S + 1)
  l_y <- numeric(S)
  delta_l <- numeric(S)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_x[1] <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- mh_sampling_step(x_curr = x[s, ], l_curr = l_x[s],
                                 l_target, ...,
                                 sampler, sampler_args = other_sampler_args,
                                 lq_sampler = lq_mh, lq_sampler_args = other_lq_args,
                                 do_checks = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_x[s+1] <- rwm_step$l_next
    y[s, ] <- rwm_step$x_prop
    l_y[s] <- rwm_step$l_prop
    delta_l[s] <- rwm_step$delta_l
    acc[s] <- rwm_step$accepted
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_x = l_x[-1],
                acc = acc, acc_rate = acc_rate,
                y = y[-1, ], l_y = l_y[-1], delta_l = delta_l[-1]))
  }
  burn_window <- 1:(burn+1)
  return(list(x = x[-burn_window, ], l_x = l_x[-burn_window],
              acc = acc[-(1:burn)], acc_rate = acc_rate,
              y = y[-(1:burn), ], l_y = l_y[-(1:burn)], delta_l = delta_l[-(1:burn)]))

}
