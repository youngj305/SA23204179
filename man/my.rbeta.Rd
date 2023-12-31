% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homework.R
\name{my.rbeta}
\alias{my.rbeta}
\title{My random beta function}
\usage{
my.rbeta(n, a, b)
}
\arguments{
\item{n}{the number of samples.}

\item{a, b}{parameters of beta distribution}
}
\value{
a random sample of size \code{n}
}
\description{
Use the acceptance-rejection method to generate a random sample from beta distribution
}
\examples{
\dontrun{
    X<-my.rbeta(1000,3,2)
}
}
