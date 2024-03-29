% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kingman_mean.R
\name{kingman_m}
\alias{kingman_m}
\title{Compute the Kingman M matrix}
\usage{
kingman_m(n = NULL, entry = NULL)
}
\arguments{
\item{n}{The number of leaves in the model}

\item{entry}{If NULL, full matrix is returned.
Else, needs to be a 2-tuple (i, j)
for the desired entry.}
}
\description{
Compute the barycentre of the F-matrices based on a formula
under the Kingman model
}
