
rwm_sampler_chain <- function(l_target, ..., scale = 1, S = 1000,
                              x_0 = NULL, x_0_u = NULL, log_lik_0 = NULL, seed = NULL,
                              custom_rw_sampler = NULL, target_names = NULL, d = NULL,
                              silent = FALSE){

  ## Preparation

  # Dimension and Scale
  stopifnot(is.numeric(scale))
  stopifnot(S>0)
  d <- d %||% ifelse(is.matrix(scale), nrow(scale), length(scale))
  if(d > 1 && !is.matrix(scale)){
    scale <- diag(scale)
    if(!silent){
      warning("Transforming scale parameter to diagonal matrix")
    }
  }


  # If the user didn't we define proposal sampler as indep. normals
  sampler <- custom_rw_sampler %||%
    ifelse(d == 1,
           function(x, scale){ rnorm(n = 1, mean = x, sd = scale) },
           function(x, scale){ mvtnorm::rmvnorm(n = 1, mean = x, sigma = scale) })

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
  log_lik_0 <- log_lik_0 %||% l_target(x_0, ...)

  # Preallocate containers
  x <- matrix(nrow = S, ncol = d, dimnames = list(NULL, target_names))
  acc <- logical(S)
  log_lik <- numeric(S)

  ## Algorithm

  # First step
  x_prop <- sampler(x_0, scale)
  log_lik_prop <- l_target(x_prop, ...)
  mstep <- mh_step(x_0, x_prop, log_lik_0, log_lik_prop, do_checks = FALSE)
  x[1, ] <- mstep$x_next
  acc[1] <- mstep$accepted
  log_lik[1] <- mstep$l_next
  if(S == 1){ # If only 1 step asked, we are done
    return(mget(c("x","acc","log_lik")))
  }

  # Resto de la cadena
  for(s in 2:S){
    x_prop <- sampler(x[s-1, ], scale)
    log_lik_prop <- l_target(x_prop, ...)
    mstep <- mh_step(x[s-1, ], x_prop, log_lik[s-1], log_lik_prop, do_checks = FALSE)
    x[s, ] <- mstep$x_next
    acc[s] <- mstep$accepted
    log_lik[s] <- mstep$l_next
  }

  return(mget(c("x","acc","log_lik")))

}
