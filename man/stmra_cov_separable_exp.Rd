% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov.R
\name{stmra_cov_separable_exp}
\alias{stmra_cov_separable_exp}
\title{Separable exp / exp spatiotemporal covariance function}
\usage{
stmra_cov_separable_exp(locs1, locs2, theta)
}
\arguments{
\item{locs1}{matrix with rows representing points and columns representing coordinates}

\item{locs2}{matrix with rows representing points and columns representing coordinates}

\item{theta}{vector of parameters (see Details)}
}
\description{
Compute pairwise covariances between two sets of locations using
a separable exponential spatiotemporal model with two exponential functions.
}
\details{
Parameter vector \code{theta} includes joint sill, spatial range, temporal range, spatial nugget, temporal nugget  (in this order).
}
