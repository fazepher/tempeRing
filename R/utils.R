
`%||%` <- function(x,y) rlang::`%||%`(x,y)

integrate_logdensity <- function(log_dens, dens_args, integrate_args = list(lower = -Inf, upper = Inf),
                                 silent = FALSE){

  dens <- function(x, log_dens, dens_args){ log_dens(match.call(dens_args)) |> exp() }

  call_args <- c(list(f = dens), log_dens = log_dens, dens_args = dens_args, integrate_args)

  result <- do.call(integrate, call_args)

  if(!silent){
    result
  }

  return(result$value)
}
