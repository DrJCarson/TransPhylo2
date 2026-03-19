// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
List num_approx_disc_multi_cpp(
    const NumericVector& grid,
    const NumericVector& off_r,
    const NumericVector& off_p,
    const NumericVector& pi,
    double w_shape,
    double w_scale,
    double ws_shape,
    double ws_scale,
    double obs_start,
    double obs_end,
    int ndemes,
    const NumericMatrix& pm
) {

  int G = grid.size();

  // Allocate outputs
  NumericMatrix omega(G, ndemes);
  NumericMatrix omega_bar(G, ndemes);
  NumericMatrix phi(G, ndemes);
  NumericMatrix pit(G, ndemes);

  // ==============================
  // Initialize first row (g = 1)
  // ==============================
  for (int d = 0; d < ndemes; d++) {
    omega(0, d)      = 1.0;
    omega_bar(0, d)  = 1.0;
    phi(0, d)        = 1.0;

    // pit[,d] = pi[d] * ( pgamma(obs_end - grid[g]) - pgamma(obs_start - grid[g]) )
    for (int g = 0; g < G; g++) {
      double x1 = obs_end  - grid[g];
      double x0 = obs_start - grid[g];
      pit(g, d) = pi[d] * ( R::pgamma(x1, ws_shape, ws_scale, 1, 0)
                              - R::pgamma(x0, ws_shape, ws_scale, 1, 0) );
    }
  }

  // =================================================
  // gamma_prob[g] = pgamma(grid[1]-grid[g]) - pgamma(grid[2]-grid[g])
  // BUT NOTE: in R, gamma_prob[1] uses grid[1] - grid[1], grid[2] - grid[1]
  // =================================================
  NumericVector gamma_prob(G, NA_REAL);

  gamma_prob[0] =
    R::pgamma(grid[0] - grid[0], w_shape, w_scale, 1, 0) -
    R::pgamma(grid[1] - grid[0], w_shape, w_scale, 1, 0);

  for (int g = 1; g < G; g++) {
    double upper = R::pgamma(grid[0] - grid[g], w_shape, w_scale, 1, 0);
    double lower = R::pgamma(grid[1] - grid[g], w_shape, w_scale, 1, 0);
    gamma_prob[g] = upper - lower;
  }

  // ft = 1 - cumsum(gamma_prob)
  NumericVector ft(G);
  double csum = 0.0;
  for (int g = 0; g < G; g++) {
    csum += gamma_prob[g];
    ft[g] = 1.0 - csum;
  }

  // Convert pm and omega[1:(g-1),] into Armadillo for matrix multiplication
  arma::mat pmA = as<arma::mat>(pm);

  // ============================================================
  // Main loop: g = 2..G (R indices 2:G, but C++ uses 1..G-1)
  // ============================================================
  for (int g = 1; g < G; g++) {

    // Extract omega[0:(g-1), ] as arma::mat   (g rows, ndemes cols)
    arma::mat omega_sub(g, ndemes);
    for (int gg = 0; gg < g; gg++) {
      for (int d = 0; d < ndemes; d++) {
        omega_sub(gg, d) = omega(gg, d);
      }
    }

    // omega_P = omega_sub %*% pm^T
    arma::mat omega_subP = omega_sub * pmA.t();        // g × ndemes

    // Loop over demes
    for (int d = 0; d < ndemes; d++) {

      // omega_bar[g,d] = ft[g] + sum( gamma_prob[g:2] * omega_P[,d] )
      // NOTE: In R, indexing gamma_prob[g:2] means descending sequence g..2 (1-based)
      double acc = ft[g];

      // gamma_prob[g], gamma_prob[g-1], ..., gamma_prob[1]
      int idx = 0;
      for (int gg = g; gg >= 1; gg--) {
        acc += gamma_prob[gg] * omega_subP(idx, d);
        idx++;
      }

      omega_bar(g, d) = acc;

      // phi[g,d] = (off_p[d] / (1 - (1-off_p[d]) * omega_bar[g,d])) ^ off_r[d]
      double denom = 1.0 - (1.0 - off_p[d]) * omega_bar(g, d);
      double tmp = off_p[d] / denom;
      phi(g, d) = std::pow(tmp, off_r[d]);

      // omega[g,d] = (1 - pit[g,d]) * phi[g,d]
      omega(g, d) = (1.0 - pit(g, d)) * phi(g, d);
    }
  }

  return List::create(
    _["omega"]      = omega,
    _["omega_bar"]  = omega_bar,
    _["phi"]        = phi,
    _["pit"]        = pit,
    _["gamma_prob"] = gamma_prob
  );
}
