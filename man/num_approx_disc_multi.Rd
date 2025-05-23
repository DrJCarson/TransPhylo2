% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{num_approx_disc_multi}
\alias{num_approx_disc_multi}
\title{Numeric approximations (discretised) of branching process}
\usage{
num_approx_disc_multi(
  grid,
  off.r,
  off.p,
  pi,
  w.shape,
  w.scale,
  ws.shape,
  ws.scale,
  obs.start,
  obs.end,
  ndemes,
  pm
)
}
\arguments{
\item{grid}{Discrete grid over which to evaluate functions}

\item{off.r}{Shape parameter for the number of offspring}

\item{off.p}{Probability parameter for the number of offspring}

\item{pi}{Probability of host being observed}

\item{w.shape}{Shape parameter of generation time distribution}

\item{w.scale}{Scale parameter of generation time distribution}

\item{ws.shape}{Shape parameter of primary sampling time distribution}

\item{ws.scale}{Scale parameter of primary sampling time distribution}

\item{obs.start}{Start date for observations}

\item{obs.end}{End date for observations}

\item{ndemes}{Number of demes}

\item{pm}{Transmission probability matrix between demes}
}
\description{
Numeric approximations (discretised) of branching process
}
