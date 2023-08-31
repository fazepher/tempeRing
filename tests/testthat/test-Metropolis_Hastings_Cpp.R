test_that("Metropolis-Hastings step works from Rcpp", {

  # We should always accept a same-place move
  x_test <- rnorm(5)
  l_test <- mvtnorm::dmvnorm(x_test, log = TRUE)
  lq_test <- mvtnorm::dmvnorm(x_test, mean = x_test, log = TRUE)
  expect_identical(
    mh_step_cpp(x_curr = x_test, x_prop = x_test,
                l_curr = l_test, l_prop = l_test,
                lq_c2p = lq_test, lq_p2c = lq_test),
    list(x_next = x_test, l_next = l_test,
         accepted = TRUE, alpha = 1, l_ratio = 0))

  # We should always accept an "up-hill" symmetrical proposal
  # (the specific values used are actually chosen differently for stress, it should still pass)
  df_test <- rpois(1, 3) + 1
  sd_test <- rgamma(1, shape = 4)
  x_test <- rt(1, df = df_test)
  y_test <- if(x_test < 0){runif(1, x_test, 0)}else{runif(1, 0, x_test)}
  expect_identical(
    mh_step_cpp(x_curr = x_test, x_prop = y_test,
                l_curr = dt(x_test, df = df_test, log = TRUE),
                l_prop = dt(y_test, df = df_test, log = TRUE),
                lq_c2p = dnorm(y_test, mean = x_test, sd = sd_test, log = TRUE),
                lq_p2c = dnorm(x_test, mean = y_test, sd = sd_test, log = TRUE)) |>
      (\(x) x[c("x_next","l_next","accepted")])(),
    list(x_next = y_test,
         l_next = dt(y_test, df = df_test, log = TRUE),
         accepted = TRUE))

  # In the full Standard Normal scenario, a proposal away from the mode should be accepted
  # with probability "its kernel"
  x_test <- rnorm(1)
  expect_equal(mh_step_cpp(x_curr = 0, x_prop = x_test,
                           l_curr = dnorm(0, log = TRUE),
                           l_prop = dnorm(x_test, log = TRUE),
                           lq_c2p = dnorm(x_test, mean = 0, log = TRUE),
                           lq_p2c = dnorm(0, mean = x_test, log = TRUE))$alpha,
               sqrt(2*pi)*dnorm(x_test))

  # Metropolis is a specific case of Metropolis-Hastings for symmetrical densities
  x_test <- rnorm(1)
  expect_identical(mh_step_cpp(x_curr = 0, x_prop = x_test,
                               l_curr = dnorm(0, log = TRUE),
                               l_prop = dnorm(x_test, log = TRUE),
                               lq_c2p = dnorm(x_test, mean = 0, log = TRUE),
                               lq_p2c = dnorm(0, mean = x_test, log = TRUE)) |>
                     (\(x) x[c("alpha","l_ratio")])(),
                   metropolis_step_cpp(x_curr = 0, x_prop = x_test,
                                       l_curr = dnorm(0, log = TRUE),
                                       l_prop = dnorm(x_test, log = TRUE)) |>
                     (\(x) x[c("alpha","l_ratio")])())

})

test_that("Metropolis-Hastings sampling step works for basic distributions", {

  # Shapiro-Wilks for basic Normals
  expect_true({
    n_test <- 1000
    mu_test <- rt(1, df = 4, ncp = 8427)
    sd_test <- sqrt(rgamma(1, shape = 997))
    x_test <- rnorm(n_test, mean = mu_test, sd = sd_test)
    l_test <- dnorm(x_test, mean = mu_test, sd = sd_test, log = TRUE)
    x_next_test <- numeric(n_test)
    for(s in 1:n_test){
      x_next_test[s] <- mh_sampling_step_cpp(x_test[s], l_test[s], l_target = dnorm,
                                             mean = mu_test, sd = sd_test, log = TRUE,
                                             sampler = function(x_c){
                                               rnorm(1, mean = x_c, sd = 2.38*sd_test)})$x_next
    }

    shapiro.test(x_next_test)$p.value > 0.025
  })

  # Kolmogorov-Smirnov for Beta-Binomial toy example
  expect_true({
    n_test <- 1000
    alpha_test <- 1.5
    beta_test <- 3.2
    p_test <- rbeta(n_test, alpha_test, beta_test)
    y_test <- sapply(p_test,function(p) rbinom(1, 5, p))
    p_post_test <- rbeta(n_test, y_test + alpha_test, 5 - y_test + beta_test)
    l_target_test <- function(p,y){
      if(p <= 0 || p >= 1){
        return(-Inf)
      }else{
        return(dbeta(p, alpha_test, beta_test, log = TRUE) + dbinom(y, 5, p, log = TRUE))
      }
    }
    p_next_test <- numeric(n_test)
    for(s in 1:n_test){
      l_curr <- l_target_test(p_test[s], y = y_test[s])
      p_curr_test <- metropolis_sampling_step_cpp(p_test[s], l_curr,
                                                  l_target_test, y = y_test[s],
                                                  sampler = function(p_c){
                                                    rnorm(1, mean = p_c, sd = 0.25)
                                                  })$x_next
      l_curr <- l_target_test(p_curr_test, y = y_test[s])
      p_next_test[s] <- metropolis_sampling_step_cpp(p_curr_test, l_curr,
                                                     l_target_test, y = y_test[s],
                                                     sampler = function(p_c){
                                                       rnorm(1, mean = p_c, sd = 0.25)
                                                     })$x_next
    }

    ks.test(p_post_test, p_next_test)$p.value > 0.025
  })


})

test_that("Random-Walk Metropolis works",{

  expect_true({
    d_test <- 5
    n_test <- 50
    logit_x_data <- matrix(rnorm(d_test*n_test), nrow = n_test, ncol = d_test)
    ltarget_logistic_reg <- function(beta, y, x = logit_x_data){
      l_prior <- sum(dnorm(beta, log = TRUE))
      eta <- as.vector(x %*% beta) # Linear predictors / Log-odds vector
      phi <- plogis(eta, log.p = TRUE) # Log-probabilities of success
      l_vero <- sum(phi + eta*(y - 1)) # Log-likelihood
      return(l_prior + l_vero)
    }

    s_test <- 500
    beta_dgp <- matrix(rnorm(d_test*s_test), nrow = d_test, ncol = s_test)
    beta_test <- matrix(rnorm(d_test*s_test), nrow = d_test, ncol = s_test)
    eta_dgp <- logit_x_data %*% beta_dgp
    eta_test <- logit_x_data %*% beta_test
    h_eta_dgp <- apply(eta_dgp, 2, function(eta_rep) c(quantile(eta_rep,c(0.05,0.5,0.95)),
                                                       sd(eta_rep),
                                                       sd(eta_rep)/mean(eta_rep)))
    h_eta_test <- apply(eta_test, 2, function(eta_rep) c(quantile(eta_rep,c(0.05,0.5,0.95)),
                                                         sd(eta_rep),
                                                         sd(eta_rep)/mean(eta_rep)))
    y_dgp <- apply(eta_dgp, 2,
                   function(eta_rep) sapply(eta_rep, function(eta) rbinom(1,1,plogis(eta))))
    y_test <- apply(eta_test, 2,
                    function(eta_rep) sapply(eta_rep, function(eta) rbinom(1,1,plogis(eta))))
    h_y_dgp <- apply(y_dgp, 2, function(y_rep) c(mean(y_rep), sd(y_rep), sd(y_rep)/mean(y_rep)))
    h_y_test <- apply(y_test, 2, function(y_rep) c(mean(y_rep), sd(y_rep), sd(y_rep)/mean(y_rep)))
    sample_dgp <- rbind(beta_dgp, h_eta_dgp, h_y_dgp)
    sample_test <- rbind(beta_test, h_eta_test, h_y_test)
    for(s in 1:s_test){
      sample_test[1:d_test, s] <- rwm_global_scale_sampler_leaner_chain(
        ltarget_logistic_reg, d = d_test, y = y_test[,s],
        global_scale = 0.35, S = 100, burn = 99,
        x_0 = beta_test[,s], silent = TRUE)$x
    }
    p_values_test <- numeric(nrow(sample_dgp))
    for(k in seq_along(p_values_test)){
      p_values_test[k] <- ks.test(sample_dgp[k, ], sample_test[k, ], exact = TRUE)$p.value
    }
    all(p_values_test > 0.025/length(p_values_test))
  })



})

test_that("Metropolis-Hastings works",{

  expect_true({
    d_test <- 5
    n_test <- 50
    logit_x_data <- matrix(rnorm(d_test*n_test), nrow = n_test, ncol = d_test)
    ltarget_logistic_reg <- function(beta, y, x = logit_x_data){
      l_prior <- sum(dnorm(beta, log = TRUE))
      eta <- as.vector(x %*% beta) # Linear predictors / Log-odds vector
      phi <- plogis(eta, log.p = TRUE) # Log-probabilities of success
      l_vero <- sum(phi + eta*(y - 1)) # Log-likelihood
      return(l_prior + l_vero)
    }
    gradient_logistic_reg <- function(beta, y, x = logit_x_data){
      eta <- as.vector(x %*% beta)
      phi <- plogis(eta, log.p = TRUE)
      grad_ll <- t(x) %*% (y - phi)
    }

    s_test <- 500
    beta_dgp <- matrix(rnorm(d_test*s_test), nrow = d_test, ncol = s_test)
    beta_test <- matrix(rnorm(d_test*s_test), nrow = d_test, ncol = s_test)
    eta_dgp <- logit_x_data %*% beta_dgp
    eta_test <- logit_x_data %*% beta_test
    h_eta_dgp <- apply(eta_dgp, 2, function(eta_rep) c(quantile(eta_rep,c(0.05,0.5,0.95)),
                                                       sd(eta_rep),
                                                       sd(eta_rep)/mean(eta_rep)))
    h_eta_test <- apply(eta_test, 2, function(eta_rep) c(quantile(eta_rep,c(0.05,0.5,0.95)),
                                                         sd(eta_rep),
                                                         sd(eta_rep)/mean(eta_rep)))
    y_dgp <- apply(eta_dgp, 2,
                   function(eta_rep) sapply(eta_rep, function(eta) rbinom(1,1,plogis(eta))))
    y_test <- apply(eta_test, 2,
                    function(eta_rep) sapply(eta_rep, function(eta) rbinom(1,1,plogis(eta))))
    h_y_dgp <- apply(y_dgp, 2, function(y_rep) c(mean(y_rep), sd(y_rep), sd(y_rep)/mean(y_rep)))
    h_y_test <- apply(y_test, 2, function(y_rep) c(mean(y_rep), sd(y_rep), sd(y_rep)/mean(y_rep)))
    sample_dgp <- rbind(beta_dgp, h_eta_dgp, h_y_dgp)
    sample_test <- rbind(beta_test, h_eta_test, h_y_test)
    for(s in 1:s_test){
      sample_test[1:d_test, s] <- mh_sampler_leaner_chain(
        ltarget_logistic_reg, d = d_test, y = y_test[,s],
        global_scale = 0.35, S = 100, burn = 99,
        x_0 = beta_test[,s], silent = TRUE)$x
    }
    p_values_test <- numeric(nrow(sample_dgp))
    for(k in seq_along(p_values_test)){
      p_values_test[k] <- ks.test(sample_dgp[k, ], sample_test[k, ], exact = TRUE)$p.value
    }
    all(p_values_test > 0.025/length(p_values_test))
  })



})
