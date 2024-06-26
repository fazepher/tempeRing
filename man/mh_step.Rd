% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Metropolis_Hastings_Step.R, R/RcppExports.R
\name{mh_step}
\alias{mh_step}
\alias{metropolis_step}
\alias{mh_sampling_step}
\alias{metropolis_sampling_step}
\alias{mh_sampling_step_cpp}
\alias{metropolis_sampling_step_cpp}
\alias{mh_step_cpp}
\title{Metropolis-Hastings Step}
\usage{
mh_step(
  x_curr,
  x_prop,
  l_curr,
  l_prop,
  lq_c2p = 0,
  lq_p2c = 0,
  do_checks = TRUE,
  full_return = TRUE
)

metropolis_step(
  x_curr,
  x_prop,
  l_curr,
  l_prop,
  do_checks = TRUE,
  full_return = TRUE
)

mh_sampling_step(
  x_curr,
  l_curr,
  l_target,
  ...,
  sampler,
  sampler_args = NULL,
  lq_sampler = NULL,
  lq_sampler_args = NULL,
  do_checks = TRUE,
  full_return = TRUE
)

metropolis_sampling_step(
  x_curr,
  l_curr,
  l_target,
  ...,
  sampler,
  sampler_args,
  do_checks = TRUE,
  full_return = TRUE
)

mh_sampling_step_cpp(
  x_curr,
  l_curr,
  l_target,
  ...,
  sampler,
  sampler_args = NULL,
  lq_sampler = NULL,
  lq_sampler_args = NULL
)

metropolis_sampling_step_cpp(
  x_curr,
  l_curr,
  l_target,
  ...,
  sampler,
  sampler_args = NULL
)

mh_step_cpp(x_curr, x_prop, l_curr, l_prop, lq_c2p = 0, lq_p2c = 0)
}
\arguments{
\item{x_curr, x_prop}{Two vectors or matrices of current and proposed states. See \code{Dimension} for
how to specify these to correctly identify the dimension of the state space.}

\item{l_curr, l_prop}{Vectors of values of log-densities, possibly up to proportionality, of the
current and proposed states.}

\item{lq_c2p, lq_p2c}{(Optional) Vectors of transition log-densities, possibly up to proportionality.
The suffix \verb{_c2p} corresponds to transitions from the \strong{c}urrent towards the \strong{p}roposed state(s),
while the suffix \verb{_p2c} contains the reversed transitions. If these values are ommited (the default),
the Hastings ratio becomes standard Metropolis.}

\item{do_checks}{If \code{TRUE} (the default), runs some preliminary checks for arguments validity.
It may be set to \code{FALSE} if one is sure that the arguments are correct and doesn't want to incurr in
slight overhead.}

\item{full_return}{If \code{TRUE} (the default), the call returns its full set of results. See  \code{Value}.}

\item{l_target}{For \code{mh_sampling_step()}, a function that computes the log-density of the target,
possibly up to proportionality. It must accept \code{x_curr} as its first argument and its returned
value must be a valid corresponding \code{l_prop}.}

\item{...}{Arguments to be passed on to the \code{l_target} function.}

\item{sampler}{For \code{mh_sampling_step()}, a sampling function whose first argument is compatible with
\code{x_curr} and returns appropriate \code{x_prop} states.}

\item{sampler_args}{A list of further arguments to be passed on to the \code{sampler}}

\item{lq_sampler}{(Optional) For \code{mh_sampling_step()}, a function that computes the transition log-densities
for the \code{sampler()}. It first two arguments must be conformable with \code{x_curr} and \code{x_prop} as
it will be used to compute the relevant \code{lq_c2p} and \code{lq_p2c} factors of the MH step.
If ommited, it is assumed that it is symmetric and, hence, it's not needed.}

\item{lq_sampler_args}{A list of further arguments to be passed on to \code{lq_sampler}.}
}
\value{
A list containing the results of the Metropolis-Hastings step:
\itemize{
\item \code{x_next}:
The next values of the chain. This object has the same structure as the input \code{x_curr}.
For each proposal, it contains the corresponding \code{x_prop} or \code{x_curr} values depending on
their acceptance or rejection.
\item \code{l_next}:
A vector of log-density values corresponding to \code{x_next} elements.
These may be re-used in the next iteration of a chain to avoid wasting resources in what
usually is an expensive computation.
\item \code{accepted}:
A vector specifying whether or not each of the proposals was accepted or rejected.
May be useful for acceptance rate monitoring.
If \code{full_return} is \code{TRUE}, then the resulting list also contains the following:
\item \code{delta_l}:
The MH log-ratio
\item The called parameters \code{x_curr}, \code{x_prop}, \code{l_curr},\code{l_prop},\code{lq_c2p},\code{lq_p2c}.
}
}
\description{
Perform a Metropolis-Hastings step on one or more states.
\itemize{
\item If we already have proposed states and all log-density factors we can call \code{mh_step()}.
\item If we wish to first sample a proposal and then give the step, we can use \code{mh_sampling_step()}.
\item Both \code{metropolis_step()} and \code{metropolis_sampling_step()} are just wrappers for symmetric
transition kernels.
}
}
\section{Rcpp}{
The functions with the suffix \verb{_cpp} internally call functions in C++ via \code{Rcpp} and should be
more efficient.
}

\section{Dimension}{
We want to perform MH steps on \code{C} states, so all log-density inputs are vectors of length \code{C}.
The only exception is for the transition log-densities that by default are set both to \code{0}
to signal that we are dealing with symmetrical Metropolis steps.
Similarly, the \code{accepted} and \code{l_next} elements of the output list are always vectors of length \code{C}.

States, on the other hand, live in a space of dimension \code{d}. There are thus 3 big combinations:
\itemize{
\item \code{C}=1: We want to perform a step on a single d-dimensional state. Without loss of generality,
\code{x_curr} and \code{x_prop} are two vectors of length \code{d}.
\item \code{C}>1, \code{d}>1: We want to perform a step on each of several multidimensional states.
Here, \code{x_curr} and \code{x_prop} are two \code{C} by \code{d} matrices, so that each state is a row vector.
\item \code{C}>1, \code{d}=1: We want to perform a step on each of several 1-dimensional states.
Here, all inputs can be vectors of length \code{C} or, alternatively,
\code{x_curr} and \code{x_prop} may be two \code{C} by \code{1} matrices.
}

The output element \code{x_next} is always of the same structure as the input \code{x_curr}.
}

\section{The Metropolis-Hastings Ratio}{
The usual form of the MH acceptance probability \eqn{\alpha = \min \lbrace 1, r \rbrace },
relies on the ratio
\deqn{ r = \pi(x_1) q(x_0|x_1) / \pi(x_0) q(x_1|x_0) }
to satisfy detail balance.
For numerical reasons, we wish to work on the log scale and the ratio becomes
\deqn{ \log(r) =  l(x_1) + lq(x_0|x_1) - l(x_0) - lq(x_1|x_0) }
Whenever the transition kernel is symmetrical (i.e. \eqn{q(x_0|x_1)=q(x_1|x_0)})
we can omit those terms from the calculation and recover the original Metropolis et. al ratio.
This is the default assumption of the \code{mh_step()} function.
\subsection{Note on proportionality}{

We want to also take advantage of the cancellation of normalizing constants so that we only need
the log-densities \emph{up to a constant of proportionality}. Note however that this refers
to the \strong{same} constant. While there are usually no problems when it comes to the target density,
care must be taken with general transition kernels as they may have different
underlying densities whose normalizing constants don't cancel.
}
}

