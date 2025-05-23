% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods_multiparm.R
\name{log_lik_ttree_multiparm_part}
\alias{log_lik_ttree_multiparm_part}
\title{Calculate transmission tree likelihood}
\usage{
log_lik_ttree_multiparm_part(
  ttree,
  grid,
  fn_list,
  off.r,
  off.p,
  pi,
  w.shape,
  w.scale,
  ws.shape,
  ws.scale,
  obs.start,
  obs.end,
  grid.delta,
  ndemes,
  pm,
  demes.prior,
  dyn_L,
  hosts
)
}
\arguments{
\item{ttree}{Transmission tree}

\item{grid}{Discrete grid over which to evaluate functions}

\item{fn_list}{Precalculated discrete approximations of exclusion probabilities}

\item{off.r}{Shape parameter for the number of offspring}

\item{off.p}{Probability parameter for the number of offspring}

\item{pi}{Probability of host being sampled}

\item{w.shape}{Shape parameter of generation time distribution}

\item{w.scale}{Scale parameter of generation time distribution}

\item{ws.shape}{Shape parameter of primary sampling time distribution}

\item{ws.scale}{Scale parameter of primary sampling time distribution}

\item{obs.start}{Start time of outbreak sampling}

\item{obs.end}{Stop time of outbreak sampling}

\item{grid.delta}{Discrete time step}

\item{ndemes}{Number of demes}

\item{pm}{Transmission probability matrix between demes}

\item{demes.prior}{Prior probability for the deme of the root host}

\item{dyn_L}{Existing likelihood matrix from the dynamic programming algorithm}

\item{hosts}{Hosts over which likelihood is calculated}
}
\description{
Calculate transmission tree likelihood
}
