test_that("Metropolis-Hastings step works", {

  # We should always accept a same-place move
  expect_identical(mh_step(x_curr = 0,
                           x_prop = 0,
                           l_curr = stats::dnorm(0, log = TRUE),
                           l_prop = stats::dnorm(0, log = TRUE),
                           lq_c2p = stats::dnorm(0, mean = 0, log = TRUE),
                           lq_p2c = stats::dnorm(0, mean = 0, log = TRUE)),
                   list(x_next = 0,
                        l_next = stats::dnorm(0, log = TRUE),
                        accepted = TRUE))

  expect_identical(mh_step(x_curr = matrix(0, nrow = 5),
                           x_prop = matrix(0, nrow = 5),
                           l_curr = rep(stats::dnorm(0, log = TRUE), 5),
                           l_prop = rep(stats::dnorm(0, log = TRUE), 5),
                           lq_c2p = rep(stats::dnorm(0, mean = 0, log = TRUE), 5),
                           lq_p2c = rep(stats::dnorm(0, mean = 0, log = TRUE), 5)),
                   list(x_next = matrix(0, nrow = 5),
                        l_next = rep(stats::dnorm(0, log = TRUE), 5),
                        accepted = rep(TRUE, 5)))

  # We should always accept an "up-hill" proposal
  expect_identical(mh_step(x_curr = 1,
                           x_prop = 0,
                           l_curr = stats::dnorm(1, log = TRUE),
                           l_prop = stats::dnorm(0, log = TRUE),
                           lq_c2p = stats::dnorm(0, mean = 1, log = TRUE),
                           lq_p2c = stats::dnorm(1, mean = 0, log = TRUE)),
                   list(x_next = 0,
                        l_next = stats::dnorm(0, log = TRUE),
                        accepted = TRUE))
  expect_identical(mh_step(x_curr = rep(1, 3),
                           x_prop = rep(0, 3),
                           l_curr = rep(stats::dnorm(1, log = TRUE), 3),
                           l_prop = rep(stats::dnorm(0, log = TRUE), 3),
                           lq_c2p = rep(stats::dnorm(0, mean = 1, log = TRUE), 3),
                           lq_p2c = rep(stats::dnorm(1, mean = 0, log = TRUE), 3)),
                   list(x_next = rep(0, 3),
                        l_next = rep(stats::dnorm(0, log = TRUE), 3),
                        accepted = rep(TRUE, 3)))

  skip_on_cran()

  # We should reject with VERY HIGH probability
  expect_identical(mh_step(x_curr = 0,
                           x_prop = 1000,
                           l_curr = stats::dnorm(0, log = TRUE),
                           l_prop = stats::dnorm(1000, log = TRUE)),
                   list(x_next = 0,
                        l_next = stats::dnorm(0, log = TRUE),
                        accepted = FALSE))

  # We should expect an acceptance rate of 25%
  expect_true(
    replicate(1000,
              mh_step(x_curr = 6, x_prop = 8,
                      l_curr = log(0.5), l_prop = log(0.25),
                      lq_c2p = log(0.7), lq_p2c = log(0.35))$accepted) |>
      mean() |>
      (\(x) abs(x - 0.25) <= qnorm(0.995)*sqrt(0.25*0.75/1000))())

  # We should expect an acceptance rate of 60%
  expect_true(
    replicate(1000,
              mh_step(x_curr = pi, x_prop = 47238,
                      l_curr = log(4), l_prop = log(4),
                      lq_c2p = log(1e-10), lq_p2c = log(6e-11), do_checks = FALSE)$accepted) |>
      mean() |>
      (\(x) abs(x - 0.60) <= qnorm(0.995)*sqrt(0.60*0.40/1000))())


})


test_that("Metropolis is a good wrapper of MH", {

  expect_identical(mh_step(x_curr = 0,
                           x_prop = 0,
                           l_curr = stats::dnorm(0, log = TRUE),
                           l_prop = stats::dnorm(0, log = TRUE),
                           lq_c2p = stats::dnorm(0, mean = 0, log = TRUE),
                           lq_p2c = stats::dnorm(0, mean = 0, log = TRUE)),
                   metropolis_step(x_curr = 0,
                                   x_prop = 0,
                                   l_curr = stats::dnorm(0, log = TRUE),
                                   l_prop = stats::dnorm(0, log = TRUE)))

  expect_identical(mh_step(x_curr = matrix(0, nrow = 5),
                           x_prop = matrix(0, nrow = 5),
                           l_curr = rep(stats::dnorm(0, log = TRUE), 5),
                           l_prop = rep(stats::dnorm(0, log = TRUE), 5),
                           lq_c2p = rep(stats::dnorm(0, mean = 0, log = TRUE), 5),
                           lq_p2c = rep(stats::dnorm(0, mean = 0, log = TRUE), 5)),
                   metropolis_step(x_curr = matrix(0, nrow = 5),
                                   x_prop = matrix(0, nrow = 5),
                                   l_curr = rep(stats::dnorm(0, log = TRUE), 5),
                                   l_prop = rep(stats::dnorm(0, log = TRUE), 5)))

})
