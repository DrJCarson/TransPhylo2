% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaries.R
\name{getIncidentCases}
\alias{getIncidentCases}
\title{Returns and/or plot numbers of sampled and unsampled cases over time}
\usage{
getIncidentCases(
  record,
  burnin = 0.5,
  numBins = 10,
  dateT = NA,
  show.plot = FALSE
)
}
\arguments{
\item{record}{Output from inferTTree function}

\item{burnin}{Proportion of the MCMC output to be discarded as burnin}

\item{numBins}{Number of time bins to compute and display incident cases}

\item{dateT}{Date when process stops (this can be Inf for fully resolved outbreaks)}

\item{show.plot}{Show a plot of incident cases over time with stacked bars}
}
\value{
List with four entries. Time is a vector of the time points. allCases is the average number of cases at each time in the posterior. sampledCases: average number of sampled cases. unsampCases: average number of unsampled cases.
}
\description{
Returns and/or plot numbers of sampled and unsampled cases over time
}
