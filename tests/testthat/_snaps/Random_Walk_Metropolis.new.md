# RWM works

    Code
      withr::with_seed(74392, {
        dplyr::mutate(dplyr::group_by(tidyr::unnest(tidyr::hoist(dplyr::mutate(
          test_scales, Chain = purrr::map(.x = Scale, function(s, ...)
            rwm_sampler_chain(scale = s, ..., silent = TRUE), l_target = lnorm)),
        Chain, "x", "acc_rate"), x), Scale), s = dplyr::row_number())
      })
    Output
      # A tibble: 50,000 x 5
      # Groups:   Scale [50]
         Scale     x acc_rate Chain                s
         <dbl> <dbl>    <dbl> <list>           <int>
       1 0.005  1.53    0.999 <named list [2]>     1
       2 0.005  1.52    0.999 <named list [2]>     2
       3 0.005  1.52    0.999 <named list [2]>     3
       4 0.005  1.53    0.999 <named list [2]>     4
       5 0.005  1.52    0.999 <named list [2]>     5
       6 0.005  1.52    0.999 <named list [2]>     6
       7 0.005  1.53    0.999 <named list [2]>     7
       8 0.005  1.53    0.999 <named list [2]>     8
       9 0.005  1.52    0.999 <named list [2]>     9
      10 0.005  1.53    0.999 <named list [2]>    10
      # ... with 49,990 more rows

