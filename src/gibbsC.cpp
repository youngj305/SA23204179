#include <Rcpp.h>
using namespace Rcpp;
//' @name rcpp_gibbs
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param a,b,n the parameters of the target bivariate density.
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//'  rcpp_gibbs(1000,2,3,10)
//'  }
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_gibbs(int N, int a,int b, int n) {
   NumericMatrix mat(N, 2);
   double x = 0, y = 0;
   for(int i = 0; i < N; i++) {
     x = rbinom(1, n, y)[0];
     y = rbeta(1, x+a,n-x+b)[0];
     mat(i, 0) = x;
     mat(i, 1) = y;
   }
   return(mat);
 }
