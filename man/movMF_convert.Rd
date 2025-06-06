% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/movMF_convert.R
\name{movMF_convert}
\alias{movMF_convert}
\title{Convert 'movMF' object}
\usage{
movMF_convert(movMF_object)
}
\arguments{
\item{movMF_object}{a movMF object of \eqn{K} components in \eqn{d} dimensions.}
}
\value{
a named list containing \describe{
\item{means}{a \eqn{(K \times d)} matrix of means}
\item{concentrations}{a \eqn{K} vector of concentrations}
\item{weights}{a \eqn{K} vector of weights}
}
}
\description{
Given an output from the movMF package's movMF function, 
convert them into the standard mixture parameter format.
}
