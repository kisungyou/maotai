% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/movMF_info.R
\name{movMF_info}
\alias{movMF_info}
\title{Extract meaningful information from the von Mises-Fisher mixture model}
\usage{
movMF_info(data, means, concentrations, weights)
}
\arguments{
\item{data}{an \eqn{(n\times d)} data matrix.}

\item{means}{an \eqn{(k\times d)} matrix of means.}

\item{concentrations}{a vector of length \eqn{k} of concentration parameters.}

\item{weights}{a vector of length \eqn{k} of mixing weights.}
}
\value{
a named list containing \describe{
\item{densities}{a vector of length \eqn{n} of the densities of the data points.}
\item{clustering}{a vector of length \eqn{n} of the hard clustering results.}
\item{loglkd}{the log-likelihood of the data.}
\item{AIC}{the Akaike information criterion.}
\item{BIC}{the Bayesian information criterion.}
}
}
\description{
Given a mixture of von Mises-Fisher distributions, this function computes 
several related quantities of the data on the unit hypersphere with respect 
to the specified model.
}
