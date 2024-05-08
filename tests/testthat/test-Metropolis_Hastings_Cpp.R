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
    n_test <- 500
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

    shapiro.test(x_next_test)$p.value > 0.01
  })

  # Kolmogorov-Smirnov for Beta-Binomial toy example
  expect_true({
    n_test <- 500
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

    ks.test(p_post_test, p_next_test)$p.value > 0.01
  })


})

test_that("Leaner Chains work for Logistic Regression",{

  # Logistic Regression test for different samplers
  ltarget_logistic_reg <- function(beta, y){
    l_prior <- sum(dnorm(beta, log = TRUE))
    eta <- as.vector(logit_x_data %*% beta) # Linear predictors / Log-odds vector
    phi <- plogis(eta, log.p = TRUE) # Log-probabilities of success
    l_vero <- sum(phi + eta*(y - 1)) # Log-likelihood
    return(l_prior + l_vero)
  }
  test_logistic_reg <- list()
  test_logistic_reg$genprior <- function(){
    rnorm(d_test)
  }
  test_logistic_reg$gendata <- function(theta){
    as.vector(logit_x_data %*% theta) |>
      sapply(function(eta) rbinom(1,1,plogis(eta)))
  }
  test_logistic_reg$test <- function(theta, dat){
    c(theta,mean(dat),
      ltarget_logistic_reg(theta,dat),
      as.vector(quantile(logit_x_data %*% theta,c(0.013,0.41,0.972))),
      sd(dat)/mean(dat))
  }

  # rwm_global_scale_sampler_leaner_chain
  d_test <- rpois(1,5)
  n_test <- rpois(1,50)
  logit_x_data <- matrix(rnorm(d_test*n_test), nrow = n_test, ncol = d_test)
  test_logistic_reg$stepMCMC <- function(theta, dat, thining){

    chain <- rwm_global_scale_sampler_leaner_chain(
      ltarget_logistic_reg, d = d_test, y = dat,
      global_scale = 0.35, S = thining, burn = thining-1,
      x_0 = theta, silent = TRUE)

    return(chain$x)

  }
  mcunit::expect_mcmc(test_logistic_reg, thinning = d_test*5)
  mcunit::expect_mcmc_reversible(test_logistic_reg, thinning = d_test*5)

  # mh_sampler_leaner_chain_cpp (RWM proposals)
  d_test <- rpois(1,5)
  n_test <- rpois(1,50)
  logit_x_data <- matrix(rnorm(d_test*n_test), nrow = n_test, ncol = d_test)
  qsampler_logistic_reg <- function(x_curr, scale){
    rnorm(length(x_curr), mean = x_curr, sd = scale)
  }
  lsampler_logistic_reg <- function(x_curr, x_prop, scale){
    l <- 0
    for(i in seq_along(x_curr)){
      l <- l + dnorm(x_prop[i], x_curr[i], scale, log = TRUE)
    }
    return(l)
  }
  test_logistic_reg$stepMCMC <- function(theta, dat, thining){

    chain <- mh_sampler_leaner_chain_cpp(
      ltarget_logistic_reg, d = d_test, y = dat,
      mh_sampler = qsampler_logistic_reg,
      other_sampler_args = list(scale = 2.38/sqrt(d_test)),
      lq_mh = lsampler_logistic_reg,
      other_lq_args = list(scale =  2.38/sqrt(d_test)),
      S = thining, burn = thining-1,
      x_0 = theta, silent = TRUE)

    return(chain$x)

  }
  mcunit::expect_mcmc(test_logistic_reg, thinning = d_test*5)
  mcunit::expect_mcmc_reversible(test_logistic_reg, thinning = d_test*5)

  # mh_sampler_leaner_chain_cpp (MALA)
  d_test <- rpois(1,5)
  n_test <- rpois(1,50)
  logit_x_data <- matrix(rnorm(d_test*n_test), nrow = n_test, ncol = d_test)
  grad_ltarget_logistic_reg <- function(beta, y, prior_sd){
    grad_l_prior <- -beta/(prior_sd*prior_sd) # Log-prior Gradient
    eta <- as.vector(logit_x_data %*% beta) # Linear predictors / Log-odds vector
    phi <- plogis(eta, log.p = TRUE) # Log-probabilities of success
    grad_l_vero <- as.vector(t(logit_x_data) %*% (y-phi)) # Log-likelihood Gradient
    return(grad_l_prior + grad_l_vero)
  }
  mala_sampler_logistic_reg <- function(x_curr, scale, y, prior_sd = 1){
    drift <- 0.5*scale^2*grad_ltarget_logistic_reg(x_curr, y, prior_sd)
    rnorm(length(x_curr), mean = x_curr + drift, sd = scale)
  }
  lmala_logistic_reg <- function(x_curr, x_prop, scale, y, prior_sd = 1){
    drift <- 0.5*scale^2*grad_ltarget_logistic_reg(x_curr, y, prior_sd)
    l <- 0
    for(i in seq_along(x_curr)){
      l <- l + dnorm(x_prop[i], x_curr[i] + drift[i], scale, log = TRUE)
    }
    return(l)
  }
  test_logistic_reg$stepMCMC <- function(theta, dat, thining){

    chain <- mh_sampler_leaner_chain_cpp(
      ltarget_logistic_reg, d = d_test, y = dat,
      mh_sampler = mala_sampler_logistic_reg,
      other_sampler_args = list(scale = 2.38/sqrt(d_test), y = dat),
      lq_mh = lmala_logistic_reg,
      other_lq_args = list(scale =  2.38/sqrt(d_test), y = dat),
      S = thining, burn = thining-1,
      x_0 = theta, silent = TRUE)

    return(chain$x)

  }
  mcunit::expect_mcmc(test_logistic_reg, thinning = d_test*5)
  mcunit::expect_mcmc_reversible(test_logistic_reg, thinning = d_test*5)

})
