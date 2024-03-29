% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapop_INLA.R
\name{BNPR_pop}
\alias{BNPR_pop}
\title{adaPop: Inference for effective population size modeling two populations jointly (INLA approximation)}
\usage{
BNPR_pop(
  tree1,
  tree2,
  samp_times1,
  samp_times2,
  lengthout = 100,
  prec_alpha = 0.01,
  prec_beta = 0.01,
  beta1_prec = 0.001,
  parSel = FALSE,
  preferential = FALSE,
  u.truncation = NA,
  l.truncation = NA,
  beta0_remove = FALSE,
  simplify = TRUE,
  derivative = FALSE,
  forward = TRUE
)
}
\arguments{
\item{tree1}{\code{phylo} object or list containing vectors of coalescent
times \code{coal_times}, sampling times \code{samp_times}, and number
sampled per sampling time \code{n_sampled}.}

\item{tree2}{Same as the previous parameters, now related to the second population.}

\item{samp_times1}{vector of sampling times related to tree1. It necessary in the case that tips are not at zero}

\item{samp_times2}{vector of sampling times related to tree2. It necessary in the case that tips are not at zero}

\item{lengthout}{number of grid points.}

\item{prec_alpha}{hyperparameter GMRF on Ne}

\item{prec_beta}{hyperparameter GMRF on gamma}

\item{beta1_prec}{hyperparameter GMRF on beta (pref sampling)}

\item{parSel}{Default is FALSE. If TRUE uses the parametric model, if FALSE, the adaptive one.}

\item{u.truncation}{It determines the upper bound the grid of the GMRF. Default is the minimum TMRCA between the two populations.}

\item{l.truncation}{It determines the lower bound the grid of the GMRF. Default is 0}

\item{beta0_remove}{It removes the intercept in the parameteric model if TRUE. Default is FALSE}

\item{simplify}{logical whether to fully bucket all Poisson points.}

\item{derivative}{logical whether to calculate estimates of the
log-derivative.}

\item{forward}{logical whether to use the finite difference approximations of
the log-derivative as a forward or backward derivative.}

\item{preferential.}{Default is FALSE. If TRUE includes also preferential sampling, if FALSE, preferential sampling is not included}
}
\value{
Effective pop size posterior, and posterior distributions of additional parameters.
\code{result} contains the INLA output, \code{data} contains the
information passed to INLA, \code{grid}  grid boundaries,
\code{x}  grid mid points \code{effpop} median effective population size estimates,
\code{effpop025} and \code{effpop975}  2.5th and 97.5th
posterior percentiles, \code{summary} contains a data.frame of the
estimates
If \code{parSel = FALSE}, \code{selInt} median gamma estimates,
\code{selInt025} and \code{selInt975}  2.5th and 97.5th
posterior percentiles, \code{selIntsummary} contains a data.frame of the
estimates.
If \code{parSel = TRUE}, \code{beta1} median beta estimates,
\code{beta1post} beta posterior distribution and \code{beta1summ}.
Same information are included also for the other parameter (if \code{beta0_remove = FALSE})
and \code{derivative} (if \code{derivative = TRUE})  log-derivative,
If \code{preferential = TRUE}, \code{sampInt} median beta, \code{sampIntmean} mean beta,
\code{sampInt975} and \code{sampInt025}97.5th and 2.5th posterior percentile beta
}
\description{
adaPop: Inference for effective population size modeling two populations jointly (INLA approximation)
}
