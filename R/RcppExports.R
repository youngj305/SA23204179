# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name rcpp_gibbs
#' @title A Gibbs sampler using Rcpp
#' @description A Gibbs sampler using Rcpp
#' @param N the number of samples
#' @param a,b,n the parameters of the target bivariate density.
#' @return a random sample of size \code{N}
#' @examples
#' \dontrun{
#'  rcpp_gibbs(1000,2,3,10)
#'  }
#' @export
rcpp_gibbs <- function(N, a, b, n) {
    .Call('_SA23204179_rcpp_gibbs', PACKAGE = 'SA23204179', N, a, b, n)
}

