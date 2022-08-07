
st_temp_step <- function(k_curr, x_curr, l_curr, l_target, ...,
                         beta_schedule, g_schedule = NULL, K = NULL){

  K <- K %||% length(beta_schedule)
  g_schedule <- g_schedule %||% rep(0, K)

  move_right <- stats::runif(1) <= 0.5
  if(k_curr == K && move_right){
    return(list("k_next" = k_curr, "acc" = FALSE, "l_next" = l_curr))
  }
  if(k_curr == 1 && !move_right){
    return(list("k_next" = k_curr, "acc" = FALSE, "l_next" = l_curr))
  }
  k_prop <- ifelse(move_right, k_curr + 1, k_curr - 1)

  delta_g <- g_schedule[k_prop] - g_schedule[k_curr]
  l_prop <- l_target(x_curr, beta_schedule[k_prop], ...)
  delta_H <- l_prop - l_curr + delta_g
  if(delta_H > 0 || log(stats::runif(1)) <= delta_H){
    return(list("k_next" = k_prop, "acc" = TRUE, "l_next" = l_prop))
  }
  return(list("k_next" = k_curr, "acc" = FALSE, "l_next" = l_curr))

}

attempt_swap <- function(x_1, x_2, beta_1, beta_2, l_1, l_2, l_target, ...){

  swaped_l <- c(l_target(x_1, beta = beta_2, ...),
                l_target(x_2, beta = beta_1, ...))
  delta_l <- sum(swaped_l) - (l_1 + l_2)

  if(delta_l > 0 || log(stats::runif(1)) <= delta_l){
    return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1), "l_next" = swaped_l))
  }
  return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2)))

}

naive_swap_move <- function(x_curr, beta_curr, k_curr, l_curr, l_target, ..., K = NULL, d = NULL){

  K <- K %||% length(k_curr)
  stopifnot(K >= 2)
  d <- d %||% ncol(x_curr)
  acc <- rep(NA, K)
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr

  # We select a single random pair of betas to attempt a swap
  b_1 <- sample(1:(K-1), 1)
  b_2 <- b_1 + 1
  # We get which "machines" have those beta indexes
  m_1 <- which(k_curr == b_1)
  m_2 <- which(k_curr == b_2)

  if(d == 1){
    nswap <- attempt_swap(x_curr[m_1], x_curr[m_2],
                          beta_curr[m_1], beta_curr[m_2],
                          l_curr[m_1], l_curr[m_2],
                          l_target, ...)
  }else{
    nswap <- attempt_swap(x_curr[1, m_1, ], x_curr[1, m_2, ],
                          beta_curr[m_1], beta_curr[m_2],
                          l_curr[m_1], l_curr[m_2],
                          l_target, ...)
  }


  m <- c(m_1, m_2)
  acc[m] <- nswap$acc
  if(nswap$acc){
    k_next[m] <- c(b_2, b_1)
    beta_next[m] <- nswap$beta_next
    l_next[m] <- nswap$l_next
  }

  return(mget(c("acc","k_next","beta_next","l_next")))

}

seo_swap_move <- function(x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                          K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){


  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)
  odd_indices <- odd_indices %||% seq(1, K, by = 2)
  even_indices <- even_indices %||% seq(2, K, by = 2)
  if(K %% 2 != 0){
    odd_indices <- odd_indices[-(K+1)/2]
  } else{
    even_indices <- even_indices[-K/2]
  }
  acc <- logical(K)
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr

  # Choose whether to swap odd or even indices with equal probability
  if(stats::runif(1) <= 0.5){
    b_1 <- odd_indices
  } else{
    b_1 <- even_indices
  }
  b_2 <- b_1 + 1

  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_curr == b_1[i])
    m_2 <- which(k_curr == b_2[i])
    if(d == 1){
      nswap <- attempt_swap(x_curr[m_1], x_curr[m_2],
                            beta_curr[m_1], beta_curr[m_2],
                            l_curr[m_1], l_curr[m_2],
                            l_target, ...)
    }else{
      nswap <- attempt_swap(x_curr[1, m_1, ], x_curr[1, m_2, ],
                            beta_curr[m_1], beta_curr[m_2],
                            l_curr[m_1], l_curr[m_2],
                            l_target, ...)
    }
    m <- c(m_1, m_2)
    if(nswap$acc){
      acc[m] <- TRUE
      k_next[m] <- c(b_2[i], b_1[i])
      beta_next[m] <- nswap$beta_next
      l_next[m] <- nswap$l_next
    }
  }

  return(mget(c("acc","k_next","beta_next","l_next")))

}

deo_swap_move <- function(c, x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                          K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  stopifnot(c >= 1)
  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)
  odd_indices <- odd_indices %||% seq(1, K, by = 2)
  even_indices <- even_indices %||% seq(2, K, by = 2)
  if(K %% 2 != 0){
    odd_indices <- odd_indices[-(K+1)/2]
  } else{
    even_indices <- even_indices[-K/2]
  }
  acc <- logical(K)
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr

  # Choose whether to swap odd or even indices deterministically based on c
  if(c %% 2 == 1){
    b_1 <- odd_indices
  } else{
    b_1 <- even_indices
  }
  b_2 <- b_1 + 1

  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_curr == b_1[i])
    m_2 <- which(k_curr == b_2[i])
    if(d == 1){
      nswap <- attempt_swap(x_curr[m_1], x_curr[m_2],
                            beta_curr[m_1], beta_curr[m_2],
                            l_curr[m_1], l_curr[m_2],
                            l_target, ...)
    }else{
      nswap <- attempt_swap(x_curr[1, m_1, ], x_curr[1, m_2, ],
                            beta_curr[m_1], beta_curr[m_2],
                            l_curr[m_1], l_curr[m_2],
                            l_target, ...)
    }
    m <- c(m_1, m_2)
    if(nswap$acc){
      acc[m] <- TRUE
      k_next[m] <- c(b_2[i], b_1[i])
      beta_next[m] <- nswap$beta_next
      l_next[m] <- nswap$l_next
    }
  }

  return(mget(c("acc","k_next","beta_next","l_next")))

}

temp_swap_move <- function(type = "deo", c = NULL,
                           x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                           K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  if(type == "deo"){
    return(deo_swap_move(c, x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                         K = K, odd_indices = odd_indices, even_indices = even_indices, d = d))
  }
  if(type == "seo"){
    return(seo_swap_move(x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                         K = K, odd_indices = odd_indices, even_indices = even_indices, d = d))
  }
  if(type == "naive"){
    return(naive_swap_move(x_curr, beta_curr, k_curr, l_curr, l_target, ..., K = K, d = d))
  }

  stop("Only recognized types of swap are: 'deo', 'seo' and 'naive'.")

}