
test_that("Multivariate normal works", {

  zeros_5 <- rep(0, 5)
  mean_5 <- c(1,0.5,0.6,0.75,0.8)
  beta_test <- 0.3

  expect_equal(lmvtnorm_temp(zeros_5, mu = mean_5),
               mvtnorm::dmvnorm(zeros_5, mean = mean_5, log = TRUE))
  expect_equal(lmvtnorm_temp(zeros_5, beta_test, mu = mean_5),
               mvtnorm::dmvnorm(zeros_5, mean = mean_5, sigma = diag(1/beta_test, 5), log = TRUE))

  x_2 <- c(2.5, 4.3)
  mean_2 <- c(2, 4)
  sigma_2 <- matrix(c(4,0.85,0.85,5), ncol = 2, byrow = TRUE)
  sigma_inv_2 <- solve(sigma_2)
  sigma_det_2 <- det(sigma_2)
  beta_2 <- beta_test^5
  expect_equal(lmvtnorm_temp(x_2, beta_2, mu = mean_2, sigma = sigma_2),
               mvtnorm::dmvnorm(x_2, mean = mean_2, sigma = sigma_2/beta_2, log = TRUE))
  expect_equal(lmvtnorm_temp(x_2, beta_2, mu = mean_2,
                             sigma_inv = sigma_inv_2, logdet_sigma = log(sigma_det_2)),
               mvtnorm::dmvnorm(x_2, mean = mean_2, sigma = sigma_2/beta_2, log = TRUE))

  # We are missing a Wishart test or other kind for normality
  # expect_equal(rmvtnorm_temp(n=1000, beta_2, mu = mean_2))
  # Also including LChol_sigma

})


test_that("Mixtures work", {

  x <- seq(-3,3,length.out = 100)
  # The more complex function at least recovers simple mixtures
  expect_equal(lmix(x,c(0.1,0.3,0.6),lnorm, mean = c(-5,0,2), sd = rep(1,3)),
               log(0.1*dnorm(x,mean=-5) + 0.3*dnorm(x,mean=0)+0.6*dnorm(x,mean=2)))
  # This difference in passed arguments is equivalent
  expect_equal(lmix(x, c(0.1,0.3,0.6),lnorm, mean = c(-5,0,2), sd = rep(1,3)),
               lmix(x, c(0.1,0.3,0.6),lnorm, mean = c(-5,0,2), shared_args = list(sd = 1)))


  # The tempered mixture with recursive function definition works
  correcto_z <- integrate(function(x) (0.5*dnorm(x, mean = -5) + 0.5*dnorm(x, mean = 5))^0.5,
                          lower = -Inf, upper = Inf)$value
  correcto <- ((0.5*dnorm(5, mean = -5) + 0.5*dnorm(5, mean = 5))^0.5)/correcto_z
  correcto_log <- log(correcto)

  expect_equal(lmix_temp(5, beta = 0.5, w = c(0.5, 0.5), ldens = lnorm, mean = c(-5, 5),
                         shared_args = list(sd = 1), log_z = NULL),
               correcto_log)


  expect_equal(ulmix_mvtnorm_temp(x = c(5,5), w = c(0.5,0.5), mu = list(c(-5,-5),c(5,5))),
               log(0.5*(mvtnorm::dmvnorm(c(5,5),mean = c(-5,-5)) +
                          mvtnorm::dmvnorm(c(5,5),mean = c(5,5)))))
  expect_equal(ulmix_mvtnorm_temp(x = c(5,5), beta = 0.3^7, w = c(0.5,0.5),
                                  mu = list(c(-5,-5),c(5,5)), sigma = list(diag(2,2),diag(4,2))),
               0.3^7*log(0.5*(mvtnorm::dmvnorm(c(5,5),mean = c(-5,-5), sigma = diag(2,2)) +
                                mvtnorm::dmvnorm(c(5,5),mean = c(5,5), sigma = diag(4,2)))))

})


test_that("Tawn et al's Skew Normal works",{

  expect_equal(lmix_hatskewnorm(rep(-45,5), mu = c(-45,-15,15,45)),
               (5*sn::dsn(x = -45, xi = c(-45,-15,15,45), omega = 1, alpha = 2, log = TRUE) +
                  rep(log(0.25),4)) |> logSumExp())

})
