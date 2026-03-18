// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double dyn_U_cpp(const NumericMatrix& ttree,
                 const NumericVector& grid,
                 const NumericMatrix& omega,
                 double w_shape,
                 double w_scale,
                 double obs_end,
                 double grid_delta,
                 int host1,   // 1-based
                 int deme1,   // 1-based
                 int host2,   // 1-based
                 int deme2,   // 1-based
                 const NumericMatrix& pm) {

  // Convert to 0-based
  int h1 = host1 - 1;
  int h2 = host2 - 1;
  int d1 = deme1 - 1;
  int d2 = deme2 - 1;

  // Infection times
  double inf_time1 = ttree(h1, 0);  // col1 in R
  double inf_time2 = ttree(h2, 0);

  // Time-index lookup
  double dt = obs_end + 1e-10 - inf_time2;
  int tidx2 = 1 + std::floor(dt / grid_delta);   // 1-based output

  // Convert to 0-based index
  int k = tidx2 - 1;

  // Safety bounds (R silently assumes grid long enough)
  if (k < 0) k = 0;
  if (k >= grid.size() - 1) k = grid.size() - 2;

  // Linear interpolation of omega
  double g_k   = grid[k];
  double g_k1  = grid[k + 1];
  double omega_k  = omega(k, d2);
  double omega_k1 = omega(k + 1, d2);

  double alpha = (inf_time2 - g_k) / (g_k1 - g_k);
  double omega_int2 = omega_k + alpha * (omega_k1 - omega_k);

  // Log-likelihood components
  double log1 = std::log(1.0 - omega_int2);

  double log2 = R::dgamma(inf_time2 - inf_time1,
                          w_shape,
                          w_scale,
                          /*log*/ true);

  double log3 = std::log(pm(d1, d2));

  return log1 + log2 + log3;
}

inline bool isFinite(double x) { return R_finite(x); }

//[[Rcpp::export]]
double dyn_T_cpp(const NumericMatrix& ttree,
                 const NumericMatrix& obs,
                 const NumericVector& grid,
                 const NumericMatrix& omega,
                 const NumericMatrix& omega_bar,
                 const NumericMatrix& pit,
                 const NumericVector& off_r,
                 const NumericVector& off_p,
                 const NumericVector& piV,
                 double ws_shape,
                 double ws_scale,
                 double obs_start,     // unused but kept for signature parity
                 double obs_end,
                 double grid_delta,
                 int host,             // 1-based
                 int deme) {           // 1-based

  int h = host - 1;
  int d = deme - 1;

  // Infection time
  double inf_time = ttree(h, 0);

  // Grid index for interpolation
  double dt = obs_end + 1e-10 - inf_time;
  int tidx = 1 + std::floor(dt / grid_delta);   // 1-based

  // Convert to 0-based safely
  int k = tidx - 1;
  int G = grid.size();
  if (k < 0) k = 0;
  if (k >= G - 1) k = G - 2;

  // Linear interpolation helper
  auto interp = [&](const NumericMatrix& M) {
    double gk = grid[k];
    double gk1 = grid[k+1];
    double Mk = M(k, d);
    double Mk1 = M(k+1, d);
    double alpha = (inf_time - gk) / (gk1 - gk);
    return Mk + alpha * (Mk1 - Mk);
  };

  double omega_int     = interp(omega);
  double omega_bar_int = interp(omega_bar);
  double pit_int       = interp(pit);

  // Count number of offspring (children)
  int inc_off = 0;
  int n = ttree.nrow();
  for (int j = 0; j < n; ++j) {
    double inf = ttree(j, 2);           // column 3: infector
    if (!NumericVector::is_na(inf) &&
        static_cast<int>(inf) == host)
      inc_off++;
  }

  // Negative binomial tail upper limit
  double q = R::qnbinom(1 - 1e-10,
                        /*size=*/off_r[d],
                                      /*prob=*/off_p[d],
                                      /*lower_tail=*/1,
                                      /*log_p=*/0);
                                      int sum_lim = static_cast<int>(std::floor(q));
                                      if (sum_lim < inc_off) sum_lim = inc_off;

                                      // Compute alpha_sum = sum_{k=inc_off..sum_lim} dnbinom(k)*choose(k,inc_off)*omega_bar^(k-inc_off)
                                      double alpha_sum = 0.0;

                                      for (int kidx = inc_off; kidx <= sum_lim; ++kidx) {
                                        double db = R::dnbinom(kidx,
                                                               off_r[d],
                                                                    off_p[d],
                                                                         /*log=*/false);

                                        double ch = R::choose(kidx, inc_off);   // choose(kidx, inc_off)

                                        int power = kidx - inc_off;
                                        double wpow = (power == 0) ? 1.0 : std::pow(omega_bar_int, power);

                                        alpha_sum += db * ch * wpow;
                                      }

                                      // If alpha_sum underflowed to 0, avoid log(0)
                                      if (alpha_sum <= 0.0 || !isFinite(alpha_sum)) {
                                        return -std::numeric_limits<double>::infinity();
                                      }

                                      // Observed or unobserved?
                                      double Tns = 0.0;
                                      if (ttree(h, 1) > 0) {        // column 2: #observations > 0
                                        // host is observed
                                        double obs_time = std::numeric_limits<double>::infinity();
                                        int m = obs.nrow();
                                        for (int j = 0; j < m; ++j) {
                                          if (static_cast<int>(obs(j, 1)) == host) {
                                            double t = obs(j, 0);
                                            if (t < obs_time) obs_time = t;
                                          }
                                        }

                                        double log_sampling =
                                          std::log(piV[d]) +
                                          R::dgamma(obs_time - inf_time, ws_shape, ws_scale, 1);

                                        Tns = log_sampling
                                          - std::log(1 - omega_int)
                                          + std::log(alpha_sum)
                                          + R::lgammafn(inc_off + 1.0);

                                      } else {
                                        // host unobserved
                                        Tns = std::log(1 - pit_int)
                                        - std::log(1 - omega_int)
                                        + std::log(alpha_sum)
                                        + R::lgammafn(inc_off + 1.0);
                                      }

                                      return Tns;
}

inline double negInf() { return -std::numeric_limits<double>::infinity(); }

// [[Rcpp::export]]
List log_lik_ttree_multiparm_cpp(const List& ttree_list,
                                 const NumericVector& grid,
                                 const List& fn_list,
                                 const NumericVector& off_r,
                                 const NumericVector& off_p,
                                 const NumericVector& piV,
                                 double w_shape, double w_scale,
                                 double ws_shape, double ws_scale,
                                 double obs_start, double obs_end,
                                 double grid_delta,
                                 int ndemes,
                                 const NumericMatrix& pm,
                                 const NumericVector& demes_prior) {

  // Extract ttree components
  NumericMatrix ttm   = ttree_list["ttree"];   // n x 3
  NumericMatrix obs   = ttree_list["obs"];     // m x 2
  IntegerVector demes = ttree_list["demes"];   // length n
  int n = ttm.nrow();

  // Precomputed arrays
  NumericMatrix omega     = fn_list["omega"];
  NumericMatrix omega_bar = fn_list["omega_bar"];
  NumericMatrix pit       = fn_list["pit"];

  // Allocate dyn_L matrix
  NumericMatrix dyn_L(n, ndemes);
  for (int i = 0; i < n; ++i)
    for (int d = 0; d < ndemes; ++d)
      dyn_L(i, d) = negInf();

  std::vector<double> lul(ndemes, negInf());

  // ---------------------------------------------------------
  // Build children lists once (avoid repeated scanning)
  // ---------------------------------------------------------
  std::vector<std::vector<int>> children_of(n + 1);

  for (int j = 0; j < n; ++j) {
    double parent = ttm(j, 2);  // col 3: infector (1-based or 0/NA)
    if (!NumericVector::is_na(parent)) {
      int p = (int)parent;
      if (p >= 1 && p <= n)
        children_of[p].push_back(j + 1);
    }
  }

  // ---------------------------------------------------------
  // host_order = order(ttree[,1], decreasing=TRUE)
  // ---------------------------------------------------------
  std::vector<int> host_order(n);
  for (int i = 0; i < n; ++i) host_order[i] = i + 1; // 1-based

  std::sort(host_order.begin(), host_order.end(),
            [&](int a, int b) {
              return ttm(a - 1, 0) > ttm(b - 1, 0);
            });

  // ---------------------------------------------------------
  // MAIN DP LOOP
  // ---------------------------------------------------------
  for (int idx = 0; idx < n; ++idx) {

    int i = host_order[idx];  // 1-based
    int i0 = i - 1;
    const auto& kids = children_of[i];
    bool is_leaf = kids.empty();

    int di = demes[i0];

    // -------------------------------------------------------
    // Case 1: Leaf
    // -------------------------------------------------------
    if (is_leaf) {

      if (di != NA_INTEGER && di > 0) {

        dyn_L(i0, di - 1) =
          dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                    off_r, off_p, piV,
                    ws_shape, ws_scale,
                    obs_start, obs_end,
                    grid_delta,
                    i, di);

      } else {

        for (int d = 1; d <= ndemes; ++d) {
          dyn_L(i0, d - 1) =
            dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                      off_r, off_p, piV,
                      ws_shape, ws_scale,
                      obs_start, obs_end,
                      grid_delta,
                      i, d);
        }
      }

    } else {
      // -------------------------------------------------------
      // Case 2: Not leaf (has children)
      // -------------------------------------------------------

      if (di != NA_INTEGER && di > 0) {

        double base = dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                                off_r, off_p, piV,
                                ws_shape, ws_scale,
                                obs_start, obs_end,
                                grid_delta,
                                i, di);

        dyn_L(i0, di - 1) = base;

        // Add children's contributions
        if (isFinite(base)) {

          for (int kid : kids) {
            int j0 = kid - 1;

            double m = negInf();
            for (int d2 = 1; d2 <= ndemes; ++d2) {

              double u = dyn_U_cpp(ttm, grid, omega,
                                   w_shape, w_scale,
                                   obs_end, grid_delta,
                                   i, di, kid, d2,
                                   pm);

              lul[d2 - 1] = u + dyn_L(j0, d2 - 1);
              if (lul[d2 - 1] > m)
                m = lul[d2 - 1];
            }

            if (!isFinite(m)) {
              dyn_L(i0, di - 1) = negInf();
              break;
            }

            double sumexp = 0.0;
            for (int d2 = 0; d2 < ndemes; ++d2)
              sumexp += std::exp(lul[d2] - m);

            dyn_L(i0, di - 1) += std::log(sumexp) + m;
          }
        }

      } else {

        // demes[i] unknown or <=0
        for (int d = 1; d <= ndemes; ++d) {

          double base = dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                                  off_r, off_p, piV,
                                  ws_shape, ws_scale,
                                  obs_start, obs_end,
                                  grid_delta,
                                  i, d);

          dyn_L(i0, d - 1) = base;

          if (isFinite(base)) {

            for (int kid : kids) {
              int j0 = kid - 1;

              double m = negInf();
              for (int d2 = 1; d2 <= ndemes; ++d2) {

                double u = dyn_U_cpp(ttm, grid, omega,
                                     w_shape, w_scale,
                                     obs_end, grid_delta,
                                     i, d, kid, d2,
                                     pm);

                lul[d2 - 1] = u + dyn_L(j0, d2 - 1);
                if (lul[d2 - 1] > m)
                  m = lul[d2 - 1];
              }

              if (!isFinite(m)) {
                dyn_L(i0, d - 1) = negInf();
                break;
              }

              double sumexp = 0.0;
              for (int d2 = 0; d2 < ndemes; ++d2)
                sumexp += std::exp(lul[d2] - m);

              dyn_L(i0, d - 1) += std::log(sumexp) + m;
            }
          }
        }
      }
    }
  }

  // ---------------------------------------------------------
  // ROOT host = last in ordered list
  // ---------------------------------------------------------
  int root = host_order.back();
  int r0 = root - 1;

  NumericVector lpl(ndemes);
  for (int d = 0; d < ndemes; ++d)
    lpl[d] = std::log(demes_prior[d]) + dyn_L(r0, d);

  double m = negInf();
  for (int d = 0; d < ndemes; ++d)
    if (lpl[d] > m) m = lpl[d];

  double loglik = negInf();
  if (isFinite(m)) {
    double sumexp = 0.0;
    for (int d = 0; d < ndemes; ++d)
      sumexp += std::exp(lpl[d] - m);
    loglik = std::log(sumexp) + m;
  }

  return List::create(
    _["loglik"] = loglik,
    _["dyn_L"] = dyn_L
  );
}


// [[Rcpp::export]]
List log_lik_ttree_multiparm_part_cpp(const List& ttree_list,
                                      const NumericVector& grid,
                                      const List& fn_list,
                                      const NumericVector& off_r,
                                      const NumericVector& off_p,
                                      const NumericVector& piV,
                                      double w_shape, double w_scale,
                                      double ws_shape, double ws_scale,
                                      double obs_start, double obs_end,
                                      double grid_delta,
                                      int ndemes,
                                      const NumericMatrix& pm,
                                      const NumericVector& demes_prior,
                                      NumericMatrix dyn_L,     // in/out
                                      IntegerVector hosts) {

  NumericMatrix ttm = ttree_list["ttree"];
  NumericMatrix obs = ttree_list["obs"];
  IntegerVector demes = ttree_list["demes"];
  int n = ttm.nrow();

  NumericMatrix omega     = fn_list["omega"];
  NumericMatrix omega_bar = fn_list["omega_bar"];
  NumericMatrix pit       = fn_list["pit"];

  // -------------------------------------------------------------
  // Build children lists once
  // -------------------------------------------------------------
  std::vector<std::vector<int>> children_of(n + 1);
  for (int j = 0; j < n; ++j) {
    double par = ttm(j, 2);     // col 3: infector
    if (!NumericVector::is_na(par)) {
      int p = (int)par;
      if (p >= 1 && p <= n)
        children_of[p].push_back(j + 1);
    }
  }

  // -------------------------------------------------------------
  // Work queue “todo”, mimic R logic:
  // always process the host with max(ttree[,1])
  // -------------------------------------------------------------
  std::vector<int> todo;
  todo.reserve(hosts.size());
  for (int h : hosts) todo.push_back(h);

  std::vector<double> lul(ndemes, negInf());

  while (!todo.empty()) {

    // Order decreasing by infection time
    std::sort(todo.begin(), todo.end(),
              [&](int a, int b) {
                return ttm(a - 1, 0) > ttm(b - 1, 0);
              });

    // Take the first
    int i = todo.front();
    todo.erase(todo.begin());
    int i0 = i - 1;

    // Reset dyn_L[i, ]
    for (int d = 0; d < ndemes; ++d)
      dyn_L(i0, d) = negInf();

    bool is_leaf = children_of[i].empty();
    int di = demes[i0];

    // ---------------------------------------------------------
    // Case 1: Leaf
    // ---------------------------------------------------------
    if (is_leaf) {

      if (di != NA_INTEGER && di > 0) {

        dyn_L(i0, di - 1) =
          dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                    off_r, off_p, piV,
                    ws_shape, ws_scale,
                    obs_start, obs_end,
                    grid_delta,
                    i, di);

      } else {
        for (int d = 1; d <= ndemes; ++d) {
          dyn_L(i0, d - 1) =
            dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                      off_r, off_p, piV,
                      ws_shape, ws_scale,
                      obs_start, obs_end,
                      grid_delta,
                      i, d);
        }
      }

      // ---------------------------------------------------------
      // Case 2: Not a leaf
      // ---------------------------------------------------------
    } else {

      const auto& kids = children_of[i];

      if (di != NA_INTEGER && di > 0) {

        double base = dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                                off_r, off_p, piV,
                                ws_shape, ws_scale,
                                obs_start, obs_end,
                                grid_delta,
                                i, di);

        dyn_L(i0, di - 1) = base;

        if (isFinite(base)) {

          for (int kid : kids) {
            int j0 = kid - 1;

            double m = negInf();
            for (int d2 = 1; d2 <= ndemes; ++d2) {

              double u = dyn_U_cpp(ttm, grid, omega,
                                   w_shape, w_scale,
                                   obs_end, grid_delta,
                                   i, di, kid, d2,
                                   pm);

              lul[d2 - 1] = u + dyn_L(j0, d2 - 1);
              if (lul[d2 - 1] > m) m = lul[d2 - 1];
            }

            if (!isFinite(m)) {
              dyn_L(i0, di - 1) = negInf();
              break;
            }

            double sumexp = 0.0;
            for (int d2 = 0; d2 < ndemes; ++d2)
              sumexp += std::exp(lul[d2] - m);

            dyn_L(i0, di - 1) += std::log(sumexp) + m;
          }
        }

      } else {

        // di unknown or <= 0
        for (int d = 1; d <= ndemes; ++d) {

          double base = dyn_T_cpp(ttm, obs, grid, omega, omega_bar, pit,
                                  off_r, off_p, piV,
                                  ws_shape, ws_scale,
                                  obs_start, obs_end,
                                  grid_delta,
                                  i, d);

          dyn_L(i0, d - 1) = base;

          if (isFinite(base)) {

            for (int kid : kids) {
              int j0 = kid - 1;

              double m = negInf();
              for (int d2 = 1; d2 <= ndemes; ++d2) {

                double u = dyn_U_cpp(ttm, grid, omega,
                                     w_shape, w_scale,
                                     obs_end, grid_delta,
                                     i, d, kid, d2,
                                     pm);

                lul[d2 - 1] = u + dyn_L(j0, d2 - 1);
                if (lul[d2 - 1] > m) m = lul[d2 - 1];
              }

              if (!isFinite(m)) {
                dyn_L(i0, d - 1) = negInf();
                break;
              }

              double sumexp = 0.0;
              for (int d2 = 0; d2 < ndemes; ++d2)
                sumexp += std::exp(lul[d2] - m);

              dyn_L(i0, d - 1) += std::log(sumexp) + m;
            }
          }
        }
      }
    }

    // ---------------------------------------------------------
    // Add parent into todo if needed
    // ---------------------------------------------------------
    double parent = ttm(i0, 2);
    if (!NumericVector::is_na(parent)) {
      int p = (int)parent;
      if (p > 0) {
        bool present = false;
        for (int v : todo)
          if (v == p) { present = true; break; }
        if (!present) todo.push_back(p);
      }
    }

  } // end while(todo)

  // ---------------------------------------------------------
  // Root: ttree[,3] == 0
  // ---------------------------------------------------------
  int root_host = -1;
  for (int i = 0; i < n; ++i) {
    double p = ttm(i, 2);
    if (NumericVector::is_na(p) || p == 0) {
      root_host = i + 1;
      break;
    }
  }

  if (root_host < 1)
    stop("No root host found (ttree[,3] == 0).");

  int r0 = root_host - 1;

  NumericVector lpl(ndemes);
  for (int d = 0; d < ndemes; ++d)
    lpl[d] = std::log(demes_prior[d]) + dyn_L(r0, d);

  double m = negInf();
  for (int d = 0; d < ndemes; ++d)
    if (lpl[d] > m) m = lpl[d];

  double loglik = negInf();
  if (isFinite(m)) {
    double sumexp = 0;
    for (int d = 0; d < ndemes; ++d)
      sumexp += std::exp(lpl[d] - m);
    loglik = std::log(sumexp) + m;
  }

  return List::create(
    _["loglik"] = loglik,
    _["dyn_L"]  = dyn_L
  );
}




//' C++ version of log_likelihood_coalescence_linear
//' @param infected_time double
//' @param final_time double
//' @param start_time double
//' @param kappa double
//' @param lambda double
//' @param branch_combs double
//' @param coalescence int (1 if event is a coalescent event, else 0)
//'
//' @return double log-likelihood increment
// [[Rcpp::export]]
double log_likelihood_coalescence_linear_cpp(double infected_time,
                                             double final_time,
                                             double start_time,
                                             double kappa,
                                             double lambda,
                                             double branch_combs,
                                             int coalescence) {

  double ll = 0.0;

  if (coalescence == 1) {

    if (lambda == 0.0) {

      // (-log(kappa)*kappa + branch_combs*final_time - branch_combs*start_time) / kappa
      ll = (-std::log(kappa) * kappa +
        branch_combs * final_time -
        branch_combs * start_time) / kappa;

    } else {

      double term1 = -std::log(lambda * (final_time - infected_time) + kappa);

      double log_start = std::log(lambda * (start_time - infected_time) + kappa);
      double log_final = std::log(lambda * (final_time - infected_time) + kappa);

      double term2 = -(branch_combs / lambda) * (log_start - log_final);

      ll = term1 + term2;
    }

  } else {  // no coalescence event

    if (branch_combs > 0) {

      if (lambda == 0.0) {

        ll = (branch_combs * final_time -
          branch_combs * start_time) / kappa;

      } else {

        double log_start = std::log(lambda * (start_time - infected_time) + kappa);
        double log_final = std::log(lambda * (final_time - infected_time) + kappa);

        ll = -(branch_combs / lambda) * (log_start - log_final);
      }

    } else {
      ll = 0.0;
    }
  }

  return ll;
}


// [[Rcpp::export]]
double log_lik_ptree_given_ctree_cpp(const List& ctree_list,
                                     double kappa,
                                     double lambda,
                                     Nullable<IntegerVector> hosts = R_NilValue)
{
  NumericMatrix ctree = ctree_list["ctree"];
  int L = ctree.nrow();

  // Determine host list
  std::vector<int> host_vec;

  if (hosts.isNull()) {
    // hosts = 1 : max(ctree[,4])
    int max_host = 0;
    for (int i = 0; i < L; ++i) {
      int h = static_cast<int>(ctree(i, 3));
      if (h > max_host) max_host = h;
    }
    host_vec.reserve(max_host);
    for (int h = 1; h <= max_host; ++h) host_vec.push_back(h);

  } else {
    IntegerVector H(hosts);
    host_vec.assign(H.begin(), H.end());
  }

  double log_lik = 0.0;

  // ----------------------------------------------------
  // For each host separately
  // ----------------------------------------------------
  for (int host : host_vec) {

    // get rows where ctree[,4] == host
    std::vector<int> rows;
    rows.reserve(16);

    for (int i = 0; i < L; ++i) {
      if (static_cast<int>(ctree(i,3)) == host) {
        rows.push_back(i);
      }
    }

    if (rows.empty())
      continue;

    // order rows by time decreasing
    std::sort(rows.begin(), rows.end(),
              [&](int a, int b) {
                return ctree(a,0) > ctree(b,0);
              });

    // infection time = ctree[which(ctree[,2]==max(rows))[1], 1]
    int max_row_index = -1;
    for (int r : rows)
      if (r > max_row_index) max_row_index = r;

    double inf_time = NA_REAL;

    for (int i = 0; i < L; ++i) {
      if (static_cast<int>(ctree(i,1)) == max_row_index) {
        // BUT R code checks ctree[,2] == max(rows)
        // i.e., column 2 (child1 index) equals max row id
        // not column 1
      }
    }

    // Correct version: find row where child1 == max(rows)
    inf_time = NA_REAL;
    for (int i = 0; i < L; ++i) {
      double child1 = ctree(i,1);  // column 2 in R
      if (!NumericVector::is_na(child1) &&
          static_cast<int>(child1) == (max_row_index + 1)) {
        inf_time = ctree(i,0);
        break;
      }
    }

    // If not found (very rare), fall back to earliest event time
    if (!R_finite(inf_time)) {
      inf_time = ctree(rows.back(), 0);
    }

    int lineages = 0;

    // ----------------------------------------------------
    // Sweep events: from most recent → oldest for this host
    // ----------------------------------------------------
    int rowsN = rows.size();
    for (int idx = 0; idx < rowsN; ++idx) {

      int row = rows[idx];
      double time1 = ctree(row, 0);

      // Leaf or internal?
      if (static_cast<int>(ctree(row,2)) == 0) {
        lineages++;
      } else {
        lineages--;
      }

      double time2;
      int is_coal = 0;

      if (idx < rowsN - 1) {
        int row2 = rows[idx + 1];
        time2 = ctree(row2, 0);

        if (static_cast<int>(ctree(row2, 2)) == 0)
          is_coal = 0;
        else
          is_coal = 1;
      } else {
        // last interval goes to inf_time
        time2 = inf_time;
        is_coal = 0;
      }

      // number of unordered lineage pairs = choose(lineages,2)
      double combs = 0.0;
      if (lineages >= 2)
        combs = (double)lineages * (lineages - 1) / 2.0;

      log_lik += log_likelihood_coalescence_linear_cpp(
        inf_time,
        time2,
        time1,
        kappa,
        lambda,
        combs,
        is_coal
      );
    }
  }

  return log_lik;
}
