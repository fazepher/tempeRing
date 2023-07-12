
# A single global scale real-valued parameter is assumed, dimension must be explicitly passed.
rwm_global_scale_sampler_leaner_chain <- function(l_target, d, ...,
                                                  global_scale = 2.38, S = 1000, burn = 0,
                                                  x_0 = NULL, x_0_u = 2, l_0 = NULL, seed = NULL,
                                                  custom_rw_sampler = NULL, more_sampler_args = NULL,
                                                  silent = FALSE){

  #--- Preparation -------------

  # If the user didn't, we define the proposal sampler as indep. normals
  if(is.null(custom_rw_sampler)){
    sampler <- function(x,scale){
      rmvtnorm_temp_cpp(n = 1, mu = x, global_scale = scale)
    }
  }else{
    sampler <- custom_rw_sampler
  }

  # Possibly set seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Get starting point
  if(is.null(x_0)){
    x_0 <- runif(d, min = -x_0_u, max = x_0_u)
  }

  # Preallocate containers
  x <- matrix(nrow = d, ncol = S + 1)
  acc <- logical(S)

  #--- Algorithm -------------

  # Initialize
  x[, 1] <- x_0
  l_x <- l_0 %||% l_target(x_0, ...)

  # Run iterations
  for(s in 1:S){
    rwm_step <- metropolis_sampling_step_list(x[, s], l_x, l_target, ...,
                                              sampler = sampler,
                                              sampler_args = c(list(global_scale), more_sampler_args))
    x[, s+1] <- rwm_step$x_next
    l_x <- rwm_step$l_next
    acc[s] <- rwm_step$accepted
  }

  #--- Post processing and result -------------

  acc_rate <- mean(acc)
  if(!silent){
    cat("Finished Sampling", sep = "\n")
    cat(paste("Acceptance Rate:", round(acc_rate,3)), sep = "\n")
  }

  if(burn == -1){
    return(list(x = x, l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  if(burn == 0){
    return(list(x = x[,-1], l_x_curr = l_x, acc = acc, acc_rate = acc_rate))
  }

  return(list(x = x[,-(1:(burn+1))], l_x_curr = l_x, acc = acc[-seq(1,burn)], acc_rate = acc_rate))

}
