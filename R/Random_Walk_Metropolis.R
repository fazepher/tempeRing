
rwm_sampler_chain <- function(l_target, ..., scale = 1, S = 1000, burn = 0,
                              x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                              custom_rw_sampler = NULL, more_sampler_args = NULL,
                              target_names = NULL, d = NULL, silent = FALSE){

#--- Preparation -------------

  # Define the dimension and sampler scale
  stopifnot(is.numeric(scale))
  d <- `%||%`(d,ifelse(is.matrix(scale), nrow(scale), length(scale)))
  if(d > 1 && !is.matrix(scale)){
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
    scale <- diag(scale)
  }

  # Checking for valid sample sizes
  stopifnot(S >= 1)
  stopifnot(0 <= burn && burn < S)

  # If the user didn't we define proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ stats::rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ mvtnorm::rmvnorm(n = 1, mean = x, sigma = scale) })

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    stopifnot(is.numeric(x_0_u) && x_0_u > 0)
    x_0 <- stats::runif(d, min = -x_0_u, max = x_0_u)
  }else{
    stopifnot(is.numeric(x_0) && length(x_0) == d)
  }

  # Preallocate containers
  x <- matrix(nrow = S + 1, ncol = d, dimnames = list(NULL, target_names))
  l_s <- numeric(S + 1)
  acc <- logical(S)

#--- Algorithm -------------

  # Initialize
  x[1, ] <- x_0
  l_s[1] <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step(x[s, ], l_s[s], l_target, ...,
                                         sampler = sampler,
                                         sampler_args = c(list(scale), more_sampler_args),
                                         do_checks = FALSE)
    x[s+1, ] <- rwm_step$x_next
    l_s[s+1] <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
  }

#--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == 0){
    return(list(x = x[-1, ], l_s = l_s[-1],
                acc = acc, acc_rate = acc_rate))
  }
  burn_window <- 1:(burn+1)
  return(list(x = x[-burn_window, ], l_s = l_s[-burn_window],
              acc = acc[-(1:burn)], acc_rate = acc_rate))

}
