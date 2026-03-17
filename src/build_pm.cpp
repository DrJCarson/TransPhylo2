// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix build_pm_cpp(const NumericVector& ext_rho, int ndemes) {
  NumericMatrix pm(ndemes, ndemes);
  for (int i = 0; i < ndemes; ++i) {
    double stay = ext_rho[i];
    double off  = (ndemes > 1) ? (1.0 - stay) / (ndemes - 1.0) : 0.0;
    for (int j = 0; j < ndemes; ++j) pm(i, j) = off;
    pm(i, i) = stay;
  }
  return pm;
}
