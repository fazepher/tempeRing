test_that("RWM works", {


  expect_true({
    true_mean <-  runif(1, min = -5, max = 5)
    samples <- rwm_sampler_chain(lnorm, mean = true_mean, sd = 3, scale = 15,
                      S = 750, burn = 50, silent = TRUE)$x
    abs(mean(samples) - true_mean) <= qnorm(0.995)*posterior::mcse_mean(samples)
    })

  skip_on_cran()

  test_scales <- c(0.005, seq(0.01, 0.1, by = 0.01),
                   seq(0.2, 1, by = 0.1),
                   seq(1.25, 3.5, by = 0.25),
                   seq(4, 8.5, by = 0.5),
                   seq(10, 100, by = 10)) |>
    tibble::tibble() |>
    setNames(c("Scale"))

  expect_snapshot(withr::with_seed(74392, {
      test_scales |>
        dplyr::mutate(
          Chain = purrr::map(.x = Scale, function(s, ...) rwm_sampler_chain(scale = s, ...),
                             l_target = lnorm)
        ) |>
        tidyr::hoist(Chain,"x","acc_rate") |>
        tidyr::unnest(x) |>
        dplyr::group_by(Scale) |>
        dplyr::mutate(s = dplyr::row_number())

  }))



})
