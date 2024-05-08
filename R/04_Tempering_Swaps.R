
####--- Original ---####

quanta_transformation <- function(x, beta_1, beta_2, mode){
  sqrt(beta_1/beta_2)*(x - mode) + mode
}

attempt_swap <- function(x_1, x_2, beta_1, beta_2, l_1, l_2, l_target, ...){

  swaped_l <- c(do.call(l_target,c(list(x = x_1, beta = beta_2), ...)),
                do.call(l_target,c(list(x = x_2, beta = beta_1), ...)))
  delta_l <- sum(swaped_l) - (l_1 + l_2)

  if(delta_l > 0 || log(runif(1)) <= delta_l){
    return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1), "l_prop" = swaped_l, "l_next" = swaped_l))
  }
  return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2), "l_prop" = swaped_l, "l_next" = c(l_1, l_2)))

}

swap_move <- function(type, j_deo,
                      x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                      K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)

  acc <- rep(-1, K)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr
  x_prop <- x_curr
  l_prop <- l_curr

  stopifnot(type %in% c("naive","seo","deo"))
  if(type == "naive"){
    # Choose the single pair to swap
    b_1 <- sample(1:(K-1),1)
  }else{

    # Even/Odd indices
    even_indices <- even_indices %||% seq(2, K, by = 2)
    odd_indices <- odd_indices %||% seq(1, K, by = 2)
    if(K %% 2 != 0){
      odd_indices <- odd_indices[-(K+1)/2]
    } else{
      even_indices <- even_indices[-K/2]
    }

    # SEO chooses at random, DEO chooses deterministically
    if(type == "seo"){
      if(runif(1) <= 0.5){
        b_1 <- odd_indices
      }else{
        b_1 <- even_indices
      }
    }else{
      if(j_deo %% 2 == 1){
        b_1 <- odd_indices
      }else{
        b_1 <- even_indices
      }
    }

  }
  b_2 <- b_1 + 1

  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_curr == b_1[i])
    m_2 <- which(k_curr == b_2[i])
    m <- c(m_1, m_2)
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
    b <- c(b_2[i], b_1[i])
    beta_next[m] <- nswap$beta_next
    l_prop[m] <- nswap$l_prop
    l_next[m] <- nswap$l_next
    acc[b] <- nswap$acc
    if(nswap$acc){
      k_next[m] <- b
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next","x_prop","l_prop")))

}

attempt_quanta <- function(mode_info,
                           x_1, x_2, beta_1, beta_2, l_1, l_2, l_target, ...,
                           unidimensional = TRUE){


  mod_1 <- modAssignment(x_1, beta_1, mode_info, assign_type = "euclidean")$A
  x_quanta_1 <- quanta_transformation(x_1, beta_1, beta_2, mode_info$modes[[mod_1]])
  mod_quanta_1 <- modAssignment(x_quanta_1, beta_2, mode_info, assign_type = "euclidean")$A

  mod_2 <- modAssignment(x_2, beta_2, mode_info, assign_type = "euclidean")$A
  x_quanta_2 <- quanta_transformation(x_2, beta_2, beta_1, mode_info$modes[[mod_2]])
  mod_quanta_2 <- modAssignment(x_quanta_2, beta_1, mode_info, assign_type = "euclidean")$A

  swaped_l <- c(do.call(l_target,c(list(x = x_quanta_1, beta = beta_2), ...)),
                do.call(l_target,c(list(x = x_quanta_2, beta = beta_1), ...)))

  if(unidimensional){
    x_curr <- c(x_1, x_2)
    x_prop <- c(x_quanta_1, x_quanta_2)
  }else{
    x_curr <- matrix(c(x_1, x_2), nrow = 2, byrow = TRUE)
    x_prop <- matrix(c(x_quanta_1, x_quanta_2), nrow = 2, byrow = TRUE)
  }
  rej_list <- list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
                   "l_next" = c(l_1, l_2), "x_next" = x_curr,
                   "l_prop" = swaped_l, "x_prop" = x_prop)


  # Early exit if the transformation is not reversible
  if(mod_1 != mod_quanta_1){
    return(rej_list)
  }
  if(mod_2 != mod_quanta_2){
    return(rej_list)
  }

  delta_l <- sum(swaped_l) - (l_1 + l_2)

  if(delta_l > 0 || log(runif(1)) <= delta_l){
    return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1),
                "l_next" = swaped_l, "x_next" = x_prop,
                "l_prop" = swaped_l, "x_prop" = x_prop))
  }
  return(rej_list)

}

quanta_move <- function(type, j_deo, mode_info,
                        x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                        K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)

  acc <- rep(-1, K)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr
  x_prop <- x_curr
  l_prop <- l_curr

  stopifnot(type %in% c("naive","seo","deo"))
  if(type == "naive"){
    # Choose the single pair to swap
    b_1 <- sample(1:(K-1),1)
  }else{

    # Even/Odd indices
    even_indices <- even_indices %||% seq(2, K, by = 2)
    odd_indices <- odd_indices %||% seq(1, K, by = 2)
    if(K %% 2 != 0){
      odd_indices <- odd_indices[-(K+1)/2]
    } else{
      even_indices <- even_indices[-K/2]
    }

    # SEO chooses at random, DEO chooses deterministically
    if(type == "seo"){
      if(runif(1) <= 0.5){
        b_1 <- odd_indices
      }else{
        b_1 <- even_indices
      }
    }else{
      if(j_deo %% 2 == 1){
        b_1 <- odd_indices
      }else{
        b_1 <- even_indices
      }
    }
  }
  b_2 <- b_1 + 1

  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_next == b_1[i])
    m_2 <- which(k_next == b_2[i])
    m <- c(m_1, m_2)
    if(d == 1){
      nswap <- attempt_quanta(mode_info,
                              x_curr[m_1], x_curr[m_2],
                              beta_curr[m_1], beta_curr[m_2],
                              l_curr[m_1], l_curr[m_2],
                              l_target, ...,
                              unidimensional = TRUE)
      x_prop[m] <- nswap$x_prop
      x_next[m] <- nswap$x_next
    }else{
      nswap <- attempt_quanta(mode_info,
                              x_curr[1, m_1, ], x_curr[1, m_2, ],
                              beta_curr[m_1], beta_curr[m_2],
                              l_curr[m_1], l_curr[m_2],
                              l_target, ...,
                              unidimensional = FALSE)
      x_prop[1, m, ] <- nswap$x_prop
      x_next[1, m, ] <- nswap$x_next
    }
    b <- c(b_2[i], b_1[i])
    beta_next[m] <- nswap$beta_next
    l_prop[m] <- nswap$l_prop
    l_next[m] <- nswap$l_next
    acc[b] <- nswap$acc
    if(nswap$acc){
      k_next[m] <- b
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next","x_prop","l_prop")))

}

# Main Swap functions
temp_swap_move <- function(type = "naive", j_deo = NULL, quanta = FALSE, mode_info = NULL,
                           x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                           K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  # QuanTA Swapping
  if(quanta){
    return(quanta_move(type, j_deo, mode_info,
                       x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                       K = K, odd_indices = odd_indices, even_indices = even_indices, d = d))
  }
  # Regular PT Swapping
  return(swap_move(type, j_deo,
                   x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                   K = K, odd_indices = odd_indices, even_indices = even_indices, d = d))

}

alps_swap_move <- function(type = "naive", j_deo = NULL, quanta_levels = NULL, mode_info = NULL,
                           x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                           K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL){

  #--- Early return of Regular PT swapping when quanta_levels is NULL -------------
  if(is.null(quanta_levels)){
    return(swap_move(type, j_deo,
                     x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                     K = K, odd_indices = odd_indices, even_indices = even_indices, d = d))
  }

  #--- Preparation -------------
  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)

  acc <- rep(-1, K)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr
  x_prop <- x_curr
  l_prop <- l_curr

  #--- Move Scheme -------------
  stopifnot(type %in% c("naive","seo","deo"))
  if(type == "naive"){
    # Choose the single pair to swap
    b_1 <- sample(1:(K-1),1)
  }else{

    # Even/Odd indices
    even_indices <- even_indices %||% seq(2, K, by = 2)
    odd_indices <- odd_indices %||% seq(1, K, by = 2)
    if(K %% 2 != 0){
      odd_indices <- odd_indices[-(K+1)/2]
    } else{
      even_indices <- even_indices[-K/2]
    }

    # SEO chooses at random, DEO chooses deterministically
    if(type == "seo"){
      if(runif(1) <= 0.5){
        b_1 <- odd_indices
      }else{
        b_1 <- even_indices
      }
    }else{
      if(j_deo %% 2 == 1){
        b_1 <- odd_indices
      }else{
        b_1 <- even_indices
      }
    }

  }
  b_2 <- b_1 + 1

  #--- Swaps -------------
  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_next == b_1[i])
    m_2 <- which(k_next == b_2[i])
    m <- c(m_1, m_2)

    # Decide if it's a QuanTA move or Regular PT swap
    if(b_1[i] %in% quanta_levels){
      if(d == 1){
        nswap <- attempt_quanta(mode_info,
                                x_curr[m_1], x_curr[m_2],
                                beta_curr[m_1], beta_curr[m_2],
                                l_curr[m_1], l_curr[m_2],
                                l_target, ...,
                                unidimensional = TRUE)
        x_prop[m] <- nswap$x_prop
        x_next[m] <- nswap$x_next
      }else{
        nswap <- attempt_quanta(mode_info,
                                x_curr[1, m_1, ], x_curr[1, m_2, ],
                                beta_curr[m_1], beta_curr[m_2],
                                l_curr[m_1], l_curr[m_2],
                                l_target, ...,
                                unidimensional = FALSE)
        x_prop[1, m, ] <- nswap$x_prop
        x_next[1, m, ] <- nswap$x_next
      }
    }else{
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
    }

    b <- c(b_2[i], b_1[i])
    beta_next[m] <- nswap$beta_next
    l_prop[m] <- nswap$l_prop
    l_next[m] <- nswap$l_next
    acc[b] <- nswap$acc
    if(nswap$acc){
      k_next[m] <- b
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next","x_prop","l_prop")))

}


####--- Llama a C++ interno ---####

attempt_quanta_list <- function(mode_info, x_1, x_2, beta_1, beta_2, l_1, l_2, l_target, ...,
                                pass_prev_mod_assignment = TRUE){

  n_modes <- length(mode_info$w)

  mod_1 <- modAssignment_cpp(x_1, beta_1,
                                         mode_info$l_target_modes, mode_info$modes,
                                         mode_info$L_inv, n_modes)
  x_quanta_1 <- quanta_transformation(x_1, beta_1, beta_2, mode_info$modes[[mod_1$A_beta]])
  mod_quanta_1 <- modAssignment_cpp(x_quanta_1, beta_2,
                                                mode_info$l_target_modes, mode_info$modes,
                                                mode_info$L_inv, n_modes)

  # Early exit if the transformation is not reversible
  if(mod_1$A_beta != mod_quanta_1$A_beta){
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
                "l_next" = c(l_1, l_2), "x_next" = list(x_1, x_2),
                "l_prop" = c(NA_real_, NA_real_), "x_prop" = list(x_quanta_1, NA_real_)))
  }

  mod_2 <- modAssignment_cpp(x_2, beta_2,
                                         mode_info$l_target_modes, mode_info$modes,
                                         mode_info$L_inv, n_modes)
  x_quanta_2 <- quanta_transformation(x_2, beta_2, beta_1, mode_info$modes[[mod_2$A_beta]])
  mod_quanta_2 <- modAssignment_cpp(x_quanta_2, beta_1,
                                                mode_info$l_target_modes, mode_info$modes,
                                                mode_info$L_inv, n_modes)

  # Early exit if the transformation is not reversible
  if(mod_2$A_beta != mod_quanta_2$A_beta){
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
                "l_next" = c(l_1, l_2), "x_next" = list(x_1, x_2),
                "l_prop" = c(NA_real_, NA_real_), "x_prop" = list(x_quanta_1, NA_real_)))
  }

  if(pass_prev_mod_assignment){
    swaped_l <- c(do.call(l_target,c(list(x = x_quanta_1, beta = beta_2,
                                          prev_mod_assign_maha = mod_quanta_1), ...)),
                  do.call(l_target,c(list(x = x_quanta_2, beta = beta_1,
                                          prev_mod_assign_maha = mod_quanta_2), ...)))
  }else{
    swaped_l <- c(do.call(l_target,c(list(x = x_quanta_1, beta = beta_2), ...)),
                  do.call(l_target,c(list(x = x_quanta_2, beta = beta_1), ...)))
  }

  x_quanta <- list(x_quanta_1, x_quanta_2)

  delta_l <- sum(swaped_l) - (l_1 + l_2)

  if(delta_l > 0 || log(runif(1)) <= delta_l){
    return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1),
                "l_next" = swaped_l, "x_next" = x_quanta,
                "l_prop" = swaped_l, "x_prop" = x_quanta))
  }
  return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
              "l_next" = c(l_1, l_2), "x_next" = list(x_1, x_2),
              "l_prop" = swaped_l, "x_prop" = x_quanta))

}

alps_swap_move_list <- function(type = "naive", deo_odd = TRUE, quanta_levels, mode_info = NULL,
                                x_curr, beta_curr, k_curr, l_curr, l_target, ...,
                                K = NULL, odd_indices = NULL, even_indices = NULL, d = NULL,
                                pass_mod_assignment = TRUE){

  #--- Preparation -------------
  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)

  acc <- rep(-1, K-1)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr
  x_prop <- x_curr
  l_prop <- l_curr

  #--- Move Scheme -------------
  stopifnot(type %in% c("naive","seo","deo"))
  if(type == "naive"){
    # Choose the single pair to swap
    b_1 <- sample(1:(K-1),1)
  }else{

    # Even/Odd indices
    even_indices <- even_indices %||% seq(2, K, by = 2)
    odd_indices <- odd_indices %||% seq(1, K, by = 2)
    if(K %% 2 != 0){
      odd_indices <- odd_indices[-(K+1)/2]
    } else{
      even_indices <- even_indices[-K/2]
    }

    # DEO chooses deterministically, SEO chooses at random
    odd_swap <- if(type == "deo"){ deo_odd }else{ runif(1) <= 0.5 }
    b_1 <- if(odd_swap){ odd_indices }else{ even_indices }

  }
  b_2 <- b_1 + 1

  #--- Swaps -------------
  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_next == b_1[i])
    m_2 <- which(k_next == b_2[i])
    m <- c(m_1, m_2)

    # Decide if it's a QuanTA move or Regular PT swap
    if(b_1[i] %in% quanta_levels){
      nswap <- attempt_quanta_list(mode_info,
                                   x_curr[[m_1]], x_curr[[m_2]],
                                   beta_curr[m_1], beta_curr[m_2],
                                   l_curr[m_1], l_curr[m_2],
                                   l_target, ...,
                                   pass_prev_mod_assignment = pass_mod_assignment)
      x_prop[m] <- nswap$x_prop
      x_next[m] <- nswap$x_next
    }else{
      nswap <- attempt_swap(x_curr[[m_1]], x_curr[[m_2]],
                            beta_curr[m_1], beta_curr[m_2],
                            l_curr[m_1], l_curr[m_2],
                            l_target, ...)
    }

    beta_next[m] <- nswap$beta_next
    l_prop[m] <- nswap$l_prop
    l_next[m] <- nswap$l_next
    acc[b_1[i]] <- nswap$acc
    if(nswap$acc){
      k_next[m] <- c(b_2[i], b_1[i])
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next","x_prop","l_prop")))

}

####--- Byproduct (FALLANDO) ---####

attempt_quanta_list_byprod <- function(mode_info, x_1, x_2, beta_1, beta_2, l_1, l_2,
                                       by_prod_1, by_prod_2,
                                       l_target_byprod, ...,
                                       pass_prev_mod_assignment = TRUE){

  n_modes <- length(mode_info$w)

  mod_1 <- modAssignment_cpp(x_1, beta_1,
                                         mode_info$l_target_modes, mode_info$modes,
                                         mode_info$L_inv, n_modes)
  x_quanta_1 <- quanta_transformation(x_1, beta_1, beta_2, mode_info$modes[[mod_1$A_beta]])
  mod_quanta_1 <- modAssignment_cpp(x_quanta_1, beta_2,
                                                mode_info$l_target_modes, mode_info$modes,
                                                mode_info$L_inv, n_modes)

  # Early exit if the transformation is not reversible
  if(mod_1$A_beta != mod_quanta_1$A_beta){
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
                "l_next" = c(l_1, l_2), "x_next" = list(x_1, x_2),
                "by_prod_next" = list(by_prod_1,by_prod_2),
                "l_prop" = c(NA_real_, NA_real_), "x_prop" = list(x_quanta_1, NA_real_)))
  }

  mod_2 <- modAssignment_cpp(x_2, beta_2,
                                         mode_info$l_target_modes, mode_info$modes,
                                         mode_info$L_inv, n_modes)
  x_quanta_2 <- quanta_transformation(x_2, beta_2, beta_1, mode_info$modes[[mod_2$A_beta]])
  mod_quanta_2 <- modAssignment_cpp(x_quanta_2, beta_1,
                                                mode_info$l_target_modes, mode_info$modes,
                                                mode_info$L_inv, n_modes)

  # Early exit if the transformation is not reversible
  if(mod_2$A_beta != mod_quanta_2$A_beta){
    return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
                "l_next" = c(l_1, l_2), "x_next" = list(x_1, x_2),
                "by_prod_next" = list(by_prod_1, by_prod_2),
                "l_prop" = c(NA_real_, NA_real_), "x_prop" = list(x_quanta_1, NA_real_)))
  }

  if(pass_prev_mod_assignment){
    swaped_l <- list(do.call(l_target_byprod,c(list(x = x_quanta_1, beta = beta_2,
                                                    prev_mod_assign_maha = mod_quanta_1), ...)),
                     do.call(l_target_byprod,c(list(x = x_quanta_2, beta = beta_1,
                                                    prev_mod_assign_maha = mod_quanta_2), ...)))
  }else{
    swaped_l <- list(do.call(l_target_byprod,c(list(x = x_quanta_1, beta = beta_2), ...)),
                     do.call(l_target_byprod,c(list(x = x_quanta_2, beta = beta_1), ...)))
  }

  x_quanta <- list(x_quanta_1, x_quanta_2)

  delta_l <- (swaped_l[[1]]$l_eval - swaped_l[[2]]$l_eval) - (l_1 + l_2)

  if(delta_l > 0 || log(runif(1)) <= delta_l){
    return(list("acc" = TRUE, "beta_next" = c(beta_2, beta_1),
                "l_next" = c(swaped_l[[1]]$l_eval, swaped_l[[2]]$l_eval),
                "x_next" = x_quanta,
                "by_prod_next" = list(swaped_l[[1]]$by_prod, swaped_l[[2]]$by_prod),
                "l_prop" = swaped_l, "x_prop" = x_quanta))
  }
  return(list("acc" = FALSE, "beta_next" = c(beta_1, beta_2),
              "l_next" = c(l_1, l_2),
              "x_next" = list(x_1, x_2),
              "by_prod_next" = list(by_prod_1, by_prod_2),
              "l_prop" = c(swaped_l[[1]]$l_eval, swaped_l[[2]]$l_eval), "x_prop" = x_quanta))

}


alps_swap_move_list_byprod <- function(type = "naive", deo_odd = TRUE, quanta_levels, mode_info = NULL,
                                       x_curr, beta_curr, k_curr, l_curr, by_prod_curr,
                                       l_target_byprod, ...,
                                       K = NULL, odd_indices = NULL, even_indices = NULL,
                                       d = NULL, dim_byprod = 1,
                                       pass_mod_assignment = TRUE){

  #--- Preparation -------------
  K <- K %||% length(k_curr)
  stopifnot(K >= 3)
  d <- d %||% ncol(x_curr)

  acc <- rep(-1, K-1)
  x_next <- x_curr
  k_next <- k_curr
  beta_next <- beta_curr
  l_next <- l_curr
  x_prop <- x_curr
  l_prop <- l_curr
  by_prod_next <- by_prod_curr

  #--- Move Scheme -------------
  stopifnot(type %in% c("naive","seo","deo"))
  if(type == "naive"){
    # Choose the single pair to swap
    b_1 <- sample(1:(K-1),1)
  }else{

    # Even/Odd indices
    even_indices <- even_indices %||% seq(2, K, by = 2)
    odd_indices <- odd_indices %||% seq(1, K, by = 2)
    if(K %% 2 != 0){
      odd_indices <- odd_indices[-(K+1)/2]
    } else{
      even_indices <- even_indices[-K/2]
    }

    # DEO chooses deterministically, SEO chooses at random
    odd_swap <- if(type == "deo"){ deo_odd }else{ runif(1) <= 0.5 }
    b_1 <- if(odd_swap){ odd_indices }else{ even_indices }

  }
  b_2 <- b_1 + 1

  #--- Swaps -------------
  for(i in seq_along(b_1)){
    # We get which "machines" have the beta indexes
    m_1 <- which(k_next == b_1[i])
    m_2 <- which(k_next == b_2[i])
    m <- c(m_1, m_2)

    # Decide if it's a QuanTA move or Regular PT swap
    if(b_1[i] %in% quanta_levels){
      nswap <- attempt_quanta_list_byprod(mode_info,
                                          x_curr[[m_1]], x_curr[[m_2]],
                                          beta_curr[m_1], beta_curr[m_2],
                                          l_curr[m_1], l_curr[m_2],
                                          by_prod_curr[[m_1]], by_prod_curr[[m_2]],
                                          l_target_byprod, ...,
                                          pass_prev_mod_assignment = pass_mod_assignment)
      x_prop[m] <- nswap$x_prop
      x_next[m] <- nswap$x_next
      by_prod_next[m] <- nswap$by_prod
    }else{
      nswap <- attempt_swap_byprod(x_curr[[m_1]], x_curr[[m_2]],
                                   beta_curr[m_1], beta_curr[m_2],
                                   l_curr[m_1], l_curr[m_2],
                                   by_prod_curr[[m_1]], by_prod_curr[[m_2]],
                                   l_target_byprod, ...)
    }

    beta_next[m] <- nswap$beta_next
    l_prop[m] <- nswap$l_prop
    l_next[m] <- nswap$l_next
    acc[b_1[i]] <- nswap$acc
    if(nswap$acc){
      k_next[m] <- c(b_2[i], b_1[i])
    }
  }

  return(mget(c("x_next","acc","k_next","beta_next","l_next","by_prod_next","x_prop","l_prop")))

}
