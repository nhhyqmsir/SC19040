% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Laplace.R
\name{Laplace.Metropolis}
\alias{Laplace.Metropolis}
\title{Implement a random walk Metropolis sampler for generating the standard Laplace distribution}
\usage{
Laplace.Metropolis(sigma, N)
}
\arguments{
\item{sigma}{the  shape paramaters of Laplace distribution}

\item{N}{the length of the series you wanted}
}
\value{
a series satisfy Laplace distribution
}
\description{
Implement a random walk Metropolis sampler for generating the standard Laplace distribution,For the increment, simulate from a normal distribution.
}
\examples{
\dontrun{
N<-1e4
sigma<-4
laplace1<-Laplace.Metropolis(sigma,N)
}
}
