% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homework.R
\name{perm_count5test}
\alias{perm_count5test}
\title{The Count 5 test for equal variances}
\usage{
perm_count5test(x, y)
}
\arguments{
\item{x}{data.}

\item{y}{data.}
}
\value{
logical number.
}
\description{
The Count 5 test for equal variances in Section 6.4 is based on the maximum number of extreme points. Example 6.15 shows that the Count 5 criterion is not applicable for unequal sample sizes. Implement a permutation test for equal variance based on the maximum number of extreme points that applies when sample sizes are not necessarily equal.
}
\examples{
\dontrun{
    n1 <- 20
    n2 <- 30
    mu1 <- mu2 <- 0
    sigma1<-sigma2 <- 1
    m <- 1000
    set.seed(1234)
    alphahat1 <- mean(replicate(m, expr={
      x <- rnorm(n1, mu1, sigma1)
      y <- rnorm(n2, mu2, sigma2)
      perm_count5test(x, y)
    }))
    alphahat1
}
}
