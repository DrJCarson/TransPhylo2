// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fast_mvrnorm_cpp(const NumericVector& mu,
                               const NumericMatrix& Sigma) {
  int d = mu.size();
  arma::vec mean = as<arma::vec>(mu);
  arma::mat cov  = as<arma::mat>(Sigma);
  arma::mat cholU = arma::chol(cov);                  // upper-triangular
  arma::vec z = arma::randn<arma::vec>(d);            // ~ N(0, I)
  arma::vec x = mean + cholU.t() * z;                 // L*z + mu
  return wrap(x);
}
