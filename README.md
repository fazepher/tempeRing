
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tempeRing

<!-- badges: start -->
<!-- badges: end -->

The goal of tempeRing is to provide an R package for my Summer Tempering
project at the Warwick Statistics CDT.

The package allows a user to run several MCMC methods, with a focus in
Parallel Tempering and Annealing for multimodal targets.

Available methods include:

1)  Random Walk Metropolis

2)  Simulated Tempering (ST)

3)  Parallel Tempering/Replica Exchange methods with several dynamics:

    1.  Standard naive implementation
    2.  Reversible dynamics via Stochastic Even-Odd swaps (SEO)
    3.  Non-reversible dynamics via Deterministic Even-Odd swaps (DEO)
        *à la* **Okabe (see also Syed)**
    4.  Accelarated Tempering for symmetric targets via QuanTA
        transformations *à la* **Tawn**
    5.  Weight-Preservation via Hessian Adjusted Targets (HAT) *à la*
        **Tawn and Roberts**

4)  Parallel Annealing with both reversible and non-reversible dynamics:

    1.  Annealed Leap-Point Sampler (ALPS)
    2.  Non-Reversible ALPS via DEO

Emphasis is on allowing us to see the methods in action on illustrative
toy examples and not in performance (e.g. everything is coded in R so
far) but an effort was made on making the R code as efficient as
possible. Faster implementations and extensions to Python and Julia are
on the horizon.
