# RWM works

    Code
      withr::with_seed(74392, {
        dplyr::mutate(dplyr::group_by(tidyr::unnest(tidyr::hoist(dplyr::mutate(
          test_scales, Chain = purrr::map(.x = Scale, function(s, ...)
            rwm_sampler_chain(scale = s, ...), l_target = lnorm)), Chain, "x",
        "acc_rate"), x), Scale), s = dplyr::row_number())
      })
    Output
      Finished Sampling
      Acceptance Rate: 0.997
      Finished Sampling
      Acceptance Rate: 0.991
      Finished Sampling
      Acceptance Rate: 0.989
      Finished Sampling
      Acceptance Rate: 0.994
      Finished Sampling
      Acceptance Rate: 0.991
      Finished Sampling
      Acceptance Rate: 0.983
      Finished Sampling
      Acceptance Rate: 0.984
      Finished Sampling
      Acceptance Rate: 0.986
      Finished Sampling
      Acceptance Rate: 0.961
      Finished Sampling
      Acceptance Rate: 0.974
      Finished Sampling
      Acceptance Rate: 0.97
      Finished Sampling
      Acceptance Rate: 0.94
      Finished Sampling
      Acceptance Rate: 0.902
      Finished Sampling
      Acceptance Rate: 0.885
      Finished Sampling
      Acceptance Rate: 0.85
      Finished Sampling
      Acceptance Rate: 0.835
      Finished Sampling
      Acceptance Rate: 0.779
      Finished Sampling
      Acceptance Rate: 0.774
      Finished Sampling
      Acceptance Rate: 0.743
      Finished Sampling
      Acceptance Rate: 0.696
      Finished Sampling
      Acceptance Rate: 0.637
      Finished Sampling
      Acceptance Rate: 0.603
      Finished Sampling
      Acceptance Rate: 0.556
      Finished Sampling
      Acceptance Rate: 0.485
      Finished Sampling
      Acceptance Rate: 0.47
      Finished Sampling
      Acceptance Rate: 0.448
      Finished Sampling
      Acceptance Rate: 0.395
      Finished Sampling
      Acceptance Rate: 0.382
      Finished Sampling
      Acceptance Rate: 0.308
      Finished Sampling
      Acceptance Rate: 0.339
      Finished Sampling
      Acceptance Rate: 0.293
      Finished Sampling
      Acceptance Rate: 0.248
      Finished Sampling
      Acceptance Rate: 0.258
      Finished Sampling
      Acceptance Rate: 0.213
      Finished Sampling
      Acceptance Rate: 0.197
      Finished Sampling
      Acceptance Rate: 0.193
      Finished Sampling
      Acceptance Rate: 0.162
      Finished Sampling
      Acceptance Rate: 0.173
      Finished Sampling
      Acceptance Rate: 0.155
      Finished Sampling
      Acceptance Rate: 0.131
      Finished Sampling
      Acceptance Rate: 0.109
      Finished Sampling
      Acceptance Rate: 0.062
      Finished Sampling
      Acceptance Rate: 0.048
      Finished Sampling
      Acceptance Rate: 0.023
      Finished Sampling
      Acceptance Rate: 0.027
      Finished Sampling
      Acceptance Rate: 0.028
      Finished Sampling
      Acceptance Rate: 0.017
      Finished Sampling
      Acceptance Rate: 0.018
      Finished Sampling
      Acceptance Rate: 0.016
      Finished Sampling
      Acceptance Rate: 0.01
      # A tibble: 50,000 x 5
      # Groups:   Scale [50]
         Scale     x acc_rate Chain                s
         <dbl> <dbl>    <dbl> <list>           <int>
       1 0.005  1.63    0.997 <named list [2]>     1
       2 0.005  1.62    0.997 <named list [2]>     2
       3 0.005  1.61    0.997 <named list [2]>     3
       4 0.005  1.61    0.997 <named list [2]>     4
       5 0.005  1.61    0.997 <named list [2]>     5
       6 0.005  1.61    0.997 <named list [2]>     6
       7 0.005  1.61    0.997 <named list [2]>     7
       8 0.005  1.61    0.997 <named list [2]>     8
       9 0.005  1.61    0.997 <named list [2]>     9
      10 0.005  1.62    0.997 <named list [2]>    10
      # ... with 49,990 more rows

