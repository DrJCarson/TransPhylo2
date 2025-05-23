% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{consTTree}
\alias{consTTree}
\title{Build a consensus transmission tree from a MCMC output}
\usage{
consTTree(record, burnin = 0.5)
}
\arguments{
\item{record}{Output from inferTTree function}

\item{burnin}{Proportion of the MCMC output to be discarded as burnin}
}
\value{
The consensus transmission tree
}
\description{
Build a consensus transmission tree from a MCMC output
}
