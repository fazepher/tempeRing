
quanta_transformation <- function(x, beta_1, beta_2, mode){
    sqrt(beta_1/beta_2)*(x - mode) + mode
}

attempt_quanta_1d <- function(mode_info, x_1, x_2, beta_1, beta_2, l_1, l_2, l_target, ...){

    mod_1 <- modAssignment(x_1, beta_1, mode_info, assign_type = "euclidean")$A
    x_quanta_1 <- quanta_transformation(x_1, beta_1, beta_2, mode_info$modes[[mod_1]])
    mod_quanta_1 <- modAssignment(x_quanta_1, beta_2, mode_info, assign_type = "euclidean")$A

    # Early exit if the transformation is not reversible
    if(mod_1 != mod_quanta_1){
        return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2), "x_next" = c(x_1, x_2)))
    }

    mod_2 <- modAssignment(x_2, beta_2, mode_info, assign_type = "euclidean")$A
    x_quanta_2 <- quanta_transformation(x_2, beta_2, beta_1, mode_info$modes[[mod_2]])
    mod_quanta_2 <- modAssignment(x_quanta_2, beta_1, mode_info, assign_type = "euclidean")$A

    # Early exit if the transformation is not reversible
    if(mod_2 != mod_quanta_2){
        return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2), "x_next" = c(x_1, x_2)))
    }

    swaped_l <- c(do.call(l_target,c(list(x = x_quanta_1, beta = beta_2), ...)),
                  do.call(l_target,c(list(x = x_quanta_2, beta = beta_1), ...)))
    delta_l <- sum(swaped_l) - (l_1 + l_2)

    if(delta_l > 0 || log(runif(1)) <= delta_l){
      return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1), "l_next" = swaped_l,
                  "x_next" = c(x_quanta_1, x_quanta_2)))
    }
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2),
                "x_next" = c(x_1, x_2)))

}

attempt_quanta_md <- function(mode_info, x_1, x_2, beta_1, beta_2, l_1, l_2, l_target, ...){

  mod_1 <- modAssignment(x_1, beta_1, mode_info, assign_type = "euclidean")$A
  x_quanta_1 <- quanta_transformation(x_1, beta_1, beta_2, mode_info$modes[[mod_1]])
  mod_quanta_1 <- modAssignment(x_quanta_1, beta_2, mode_info, assign_type = "euclidean")$A

  # Early exit if the transformation is not reversible
  if(mod_1 != mod_quanta_1){
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2),
                "x_next" = matrix(c(x_1, x_2), nrow = 2, byrow = TRUE)))
  }

  mod_2 <- modAssignment(x_2, beta_2, mode_info, assign_type = "euclidean")$A
  x_quanta_2 <- quanta_transformation(x_2, beta_2, beta_1, mode_info$modes[[mod_2]])
  mod_quanta_2 <- modAssignment(x_quanta_2, beta_1, mode_info, assign_type = "euclidean")$A

  # Early exit if the transformation is not reversible
  if(mod_2 != mod_quanta_2){
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2),
                "x_next" = matrix(c(x_1, x_2), nrow = 2, byrow = TRUE)))
  }

  swaped_l <- c(do.call(l_target,c(list(x = x_quanta_1, beta = beta_2), ...)),
                do.call(l_target,c(list(x = x_quanta_2, beta = beta_1), ...)))
  delta_l <- sum(swaped_l) - (l_1 + l_2)

  if(delta_l > 0 || log(runif(1)) <= delta_l){
    return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1), "l_next" = swaped_l,
                "x_next" = matrix(c(x_quanta_1, x_quanta_2), nrow = 2, byrow = TRUE)))
  }
  return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_next" = c(l_1, l_2),
              "x_next" = matrix(c(x_1, x_2), nrow = 2, byrow = TRUE)))

}


naive_quanta_move <- function(mode_info, x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                              K = NULL, d = NULL){

  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)

  acc <- rep(-1, K)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr

  # Choose with replacement which indexes we attempt to swap
  b_1 <- sample(1:(K-1),1)
  b_2 <- b_1 + 1

  # We get which "machines" have the beta indexes
  m_1 <- which(k_next == b_1)
  m_2 <- which(k_next == b_2)
  m <- c(m_1, m_2)


  if(d == 1){
    nswap <- attempt_quanta_1d(mode_info,
                               x_next[m_1], x_next[m_2],
                               beta_next[m_1], beta_next[m_2],
                               l_next[m_1], l_next[m_2],
                               l_target, ...)
    x_next[m] <- nswap$x_next
  }else{
    nswap <- attempt_quanta_md(mode_info,
                               x_next[1, m_1, ], x_next[1, m_2, ],
                               beta_next[m_1], beta_next[m_2],
                               l_next[m_1], l_next[m_2],
                               l_target, ...)
    x_next[1, m, ] <- nswap$x_next
  }
  b <- c(b_2, b_1)
  beta_next[m] <- nswap$beta_next
  l_next[m] <- nswap$l_next
  if(nswap$acc){
    k_next[m] <- b
    acc[b] <- TRUE
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next")))

}

seo_quanta_move <- function(mode_info, x_curr, beta_curr, k_curr, l_curr, l_target, ...,
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
  acc <- rep(-1,K)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr

  # Choose whether to swap odd or even indices with equal probability
  if(runif(1) <= 0.5){
    b_1 <- odd_indices
  } else{
    b_1 <- even_indices
  }
  b_2 <- b_1 + 1

  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_next == b_1[i])
    m_2 <- which(k_next == b_2[i])
    m <- c(m_1, m_2)

    if(d == 1){
      nswap <- attempt_quanta_1d(mode_info,
                                 x_next[m_1], x_next[m_2],
                                 beta_next[m_1], beta_next[m_2],
                                 l_next[m_1], l_next[m_2],
                                 l_target, ...)
      x_next[m] <- nswap$x_next
    }else{
      nswap <- attempt_quanta_md(mode_info,
                                 x_next[1, m_1, ], x_next[1, m_2, ],
                                 beta_next[m_1], beta_next[m_2],
                                 l_next[m_1], l_next[m_2],
                                 l_target, ...)
      x_next[1, m, ] <- nswap$x_next
    }
    b <- c(b_2[i], b_1[i])
    beta_next[m] <- nswap$beta_next
    l_next[m] <- nswap$l_next
    if(nswap$acc){
      k_next[m] <- b
      acc[b] <- TRUE
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next")))

}

deo_quanta_move <- function(mode_info, j_deo, x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                            K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  stopifnot(j_deo >= 1)
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
  acc <- rep(-1, K)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr

  # Choose whether to swap odd or even indices deterministically based on c
  if(j_deo %% 2 == 1){
    b_1 <- odd_indices
  } else{
    b_1 <- even_indices
  }
  b_2 <- b_1 + 1

  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_next == b_1[i])
    m_2 <- which(k_next == b_2[i])
    m <- c(m_1, m_2)

    if(d == 1){
      nswap <- attempt_quanta_1d(mode_info,
                                 x_next[m_1], x_next[m_2],
                                 beta_next[m_1], beta_next[m_2],
                                 l_next[m_1], l_next[m_2],
                                 l_target, ...)
      x_next[m] <- nswap$x_next
    }else{
      nswap <- attempt_quanta_md(mode_info,
                                 x_next[1, m_1, ], x_next[1, m_2, ],
                                 beta_next[m_1], beta_next[m_2],
                                 l_next[m_1], l_next[m_2],
                                 l_target, ...)
      x_next[1, m, ] <- nswap$x_next
    }
    b <- c(b_2[i], b_1[i])
    beta_next[m] <- nswap$beta_next
    l_next[m] <- nswap$l_next
    if(nswap$acc){
      k_next[m] <- b
      acc[b] <- TRUE
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next")))

}
