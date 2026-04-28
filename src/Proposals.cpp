// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>

#include <algorithm>  // std::sort
#include <vector>     // std::vector
#include <numeric>    // std::accumulate
#include <cmath>      // std::log, std::exp

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
List extract_host_minitree_cpp(const List& ctree_list,
                               int host)
{
  NumericMatrix ctree = ctree_list["ctree"];
  const int L = ctree.nrow();

  if (host < 1)
    stop("host must be >= 1");

  // ------------------------------------------------------------
  // rows_host: all rows belonging to host
  // ------------------------------------------------------------
  std::vector<int> rows_host;
  for (int i = 0; i < L; ++i) {
    if ((int)ctree(i,3) == host)
      rows_host.push_back(i);
  }

  if (rows_host.empty())
    stop("Host has no rows in ctree");

  const int R = rows_host.size();

  // ------------------------------------------------------------
  // leaves: any row with child2 == 0
  // (R condition: ctree[,3] == 0 & ctree[,4] == host)
  // ------------------------------------------------------------
  std::vector<int> leaves;
  for (int r : rows_host) {
    if (ctree(r,2) == 0)
      leaves.push_back(r);
  }

  if (leaves.empty())
    stop("Host has no leaves");

  const int Lf = leaves.size();

  // ------------------------------------------------------------
  // type and host2 (leaf classification)
  // type:
  //   1 = observation leaf (child1 == 0)
  //   2 = transmission leaf (child1 > 0)
  // ------------------------------------------------------------
  IntegerVector type(Lf);
  IntegerVector host2(Lf);

  for (int i = 0; i < Lf; ++i) {
    int row = leaves[i];

    // observation leaf
    if (ctree(row,1) == 0) {
      type[i]  = 1;
      host2[i] = host;

      // transmission leaf
    } else {
      type[i] = 2;

      int parent = (int)ctree(row,1) - 1;
      if (parent < 0 || parent >= L)
        stop("Invalid transmission parent index");

      host2[i] = (int)ctree(parent,3);
    }
  }

  // ------------------------------------------------------------
  // unique_hosts, hosts_count, hosts_type
  // ------------------------------------------------------------
  IntegerVector unique_hosts = sort_unique(host2);
  const int H = unique_hosts.size();

  IntegerVector hosts_count(H);
  IntegerVector hosts_type(H);

  for (int h = 0; h < H; ++h) {
    int uh = unique_hosts[h];
    int count = 0;
    int first_type = -1;

    for (int i = 0; i < Lf; ++i) {
      if (host2[i] == uh) {
        ++count;
        if (first_type < 0)
          first_type = type[i];
      }
    }

    hosts_count[h] = count;
    hosts_type[h]  = first_type;
  }

  // ------------------------------------------------------------
  // down_hosts and interval
  // ------------------------------------------------------------
  NumericMatrix down_hosts(R, H);
  NumericMatrix interval(R, 2);

  for (int l = 0; l < Lf; ++l) {

    int r = leaves[l];
    int h = host2[l];

    // find local host index
    int uh_idx = -1;
    for (int k = 0; k < H; ++k) {

      if (unique_hosts[k] == h)
        uh_idx = k;

    }

    if (uh_idx < 0)
      stop("Host lookup failed");

    // walk backward until infection event
    while (true) {

      // locate r within rows_host
      int r_idx = -1;
      for (int i = 0; i < R; ++i) {

        if (rows_host[i] == r)
          r_idx = i;

      }

      if (r_idx < 0)
        stop("Row not found in rows_host");

      // update downstream counts
      down_hosts(r_idx, uh_idx) += 1;

      // find ancestor (exactly one expected)
      int anc = -1;
      for (int i = 0; i < L; ++i) {
        if ((int)ctree(i,1) == r + 1 || (int)ctree(i,2) == r + 1) {
          anc = i;
          break;
        }
      }

      if (anc < 0)
        stop("Ancestor not found");

      interval(r_idx,0) = ctree(anc,0);
      interval(r_idx,1) = ctree(r,0);

      // stop at infection event
      if ((int)ctree(r,3) != (int)ctree(anc,3))
        break;

      r = anc;
    }
  }

  return List::create(
    _["rows_host"]    = wrap(rows_host),
    _["leaves"]       = wrap(leaves),
    _["type"]         = type,
    _["host2"]        = host2,
    _["unique_hosts"] = unique_hosts,
    _["hosts_count"]  = hosts_count,
    _["hosts_type"]   = hosts_type,
    _["down_hosts"]   = down_hosts,
    _["interval"]     = interval
  );
}


void combn_recursive(int offset, int k, int n,
                     std::vector<int>& curr,
                     std::vector< std::vector<int> >& out) {

  if (k == 0) {
    out.push_back(curr);
    return;
  }

  for (int i = offset; i <= n - k; ++i) {
    curr.push_back(i);
    combn_recursive(i + 1, k - 1, n, curr, out);
    curr.pop_back();
  }
}

//' @export
// [[Rcpp::export]]
std::vector< std::vector< Rcpp::List > > enumerate_transmission_clusters_cpp(
    const NumericMatrix& down_hosts,
    const IntegerVector& hosts_count,
    const IntegerVector& hosts_type,
    const IntegerVector& unique_hosts,
    const NumericMatrix& interval) {

  const int R = interval.nrow();     // number of branches
  const int H = down_hosts.ncol();   // number of hosts

  // ------------------------------------------------------
  // Collect unique time cutpoints
  // ------------------------------------------------------
  std::vector<double> cutpoints;
  cutpoints.reserve(2 * R);
  for (int r = 0; r < R; ++r) {
    cutpoints.push_back(interval(r, 0));
    cutpoints.push_back(interval(r, 1));
  }
  std::sort(cutpoints.begin(), cutpoints.end());
  cutpoints.erase(std::unique(cutpoints.begin(), cutpoints.end()),
                  cutpoints.end());

  // tr_lin[lb] = configurations with lb transmitted lineages
  std::vector< std::vector< List > > tr_lin(R + 1);

  // ------------------------------------------------------
  // Loop over time intervals
  // ------------------------------------------------------
  for (std::size_t ti = 0; ti + 1 < cutpoints.size(); ++ti) {

    const double start = cutpoints[ti];
    const double end   = cutpoints[ti + 1];
    if (end <= start)
      continue;

    // --------------------------------------------------
    // Branches overlapping this time interval
    // --------------------------------------------------
    std::vector<int> branches_todo;
    for (int r = 0; r < R; ++r) {
      if (!(interval(r, 1) <= start || interval(r, 0) >= end))
        branches_todo.push_back(r);
    }

    // --------------------------------------------------
    // Construct clusters
    // --------------------------------------------------
    std::vector< List > clusters;

    while (!branches_todo.empty()) {

      // start with a single branch
      std::vector<int> branches_inc;
      branches_inc.push_back(branches_todo.front());
      branches_todo.erase(branches_todo.begin());

      std::vector<int> counts_inc(H, 0);
      for (int h = 0; h < H; ++h)
        counts_inc[h] += down_hosts(branches_inc[0], h);

      bool changed = true;
      while (changed) {
        changed = false;

        for (auto it = branches_todo.begin();
             it != branches_todo.end(); ) {

          bool intersects = false;
          for (int h = 0; h < H; ++h) {
            if (counts_inc[h] > 0 &&
                down_hosts(*it, h) > 0) {
              intersects = true;
              break;
            }
          }

          if (intersects) {
            // absorb branch into cluster
            branches_inc.push_back(*it);
            for (int h = 0; h < H; ++h)
              counts_inc[h] += down_hosts(*it, h);
            it = branches_todo.erase(it);
            changed = true;
          } else {
            ++it;
          }
        }
      }

      // --------------------------------------------------
      // Determine affected hosts
      // --------------------------------------------------
      std::vector<int> hosts_inc;
      std::vector<int> obs_inc;

      for (int h = 0; h < H; ++h) {
        if (counts_inc[h] > 0) {
          hosts_inc.push_back(unique_hosts[h]);
          if (hosts_type[h] == 1)
            obs_inc.push_back(unique_hosts[h]);
        }
      }

      // At most one observed host
      if (obs_inc.size() > 1)
        continue;

      // Must include full leaf count for each host
      bool valid = true;
      for (int h = 0; h < H; ++h) {
        if (counts_inc[h] > 0 &&
            counts_inc[h] != hosts_count[h]) {
          valid = false;
          break;
        }
      }
      if (!valid)
        continue;

      // Compute bin (over host indices)
      int bin = 0;
      for (int h = 0; h < H; ++h) {
        if (counts_inc[h] > 0)
          bin |= (1 << h);

      }

      // Store cluster
      clusters.push_back(List::create(
          _["branches"] = branches_inc,   // LOCAL indices
          _["bin"]      = bin,
          _["hosts"]    = hosts_inc,
          _["obs"]      = obs_inc
      ));
    }

    // --------------------------------------------------
    // Combine clusters into transmission configurations
    // --------------------------------------------------
    const int K = clusters.size();
    if (K == 0)
      continue;

    // generate combinations
    for (int k = 0; k < K; ++k) {
      // singleton combinations first (k = 1)
      List cl = clusters[k];
      IntegerVector br = cl["branches"];
      const int lb = br.size();

      if (lb > 0) {
        tr_lin[lb].push_back(List::create(
            _["branches"] = br,
            _["bin"]      = cl["bin"],
                              _["start"]    = start,
                              _["end"]      = end,
                              _["length"]   = end - start
        ));
      }
    }

    // now combinations of size >= 2
    for (int a = 0; a < K; ++a) {
      for (int b = a + 1; b < K; ++b) {

        IntegerVector br1 = clusters[a]["branches"];
        IntegerVector br2 = clusters[b]["branches"];

        std::vector<int> merged;
        merged.reserve(br1.size() + br2.size());

        for (int x : br1) merged.push_back(x);
        for (int x : br2) merged.push_back(x);

        const int lb = merged.size();
        if (lb == 0)
          continue;

        int bin = as<int>(clusters[a]["bin"]) |
          as<int>(clusters[b]["bin"]);

        tr_lin[lb].push_back(List::create(
            _["branches"] = merged,
            _["bin"]      = bin,
            _["start"]    = start,
            _["end"]      = end,
            _["length"]   = end - start
        ));
      }
    }
  }

  // ------------------------------------------------------
  // Safety check: invariants
  // ------------------------------------------------------
  for (std::size_t lb = 0; lb < tr_lin.size(); ++lb) {
    for (const auto& entry : tr_lin[lb]) {

      IntegerVector branches = entry["branches"];

      if (branches.size() != static_cast<R_xlen_t>(lb)) {
        stop("Invariant violated: branches.size() != lb");
      }

      for (int b : branches) {
        if (b < 0 || b >= R) {
          stop("Invariant violated: branch index out of bounds");
        }
      }
    }
  }

  return tr_lin;

}


int sample_index(const Rcpp::NumericVector& prob) {
  double u = R::runif(0.0, 1.0);
  double cum = 0.0;

  for (int i = 0; i < prob.size(); ++i) {
    cum += prob[i];
    if (u <= cum)
      return i + 1;  // 1-based to match R
  }

  // Fallback for numeric roundoff
  return prob.size();
}


//' @export
// [[Rcpp::export]]
Rcpp::List sample_transmission_cpp(const Rcpp::List& tr_lin)
{
  using namespace Rcpp;

  const int L = tr_lin.size();

  // ------------------------------------------------------------
  // Compute total length for each lineage count
  // ------------------------------------------------------------
  NumericVector len_lin(L, 0.0);

  for (int l = 0; l < L; ++l) {
    List configs = tr_lin[l];
    double total = 0.0;

    for (int j = 0; j < configs.size(); ++j)
      total += as<double>(as<List>(configs[j])["length"]);

    len_lin[l] = total;
  }

  // Feasible lineage counts

  std::vector<int> poss_lin;
  for (int i = 0; i < len_lin.size(); ++i) {
    if (len_lin[i] > 0)
      poss_lin.push_back(i);  // 0-based
  }

  if (poss_lin.empty())
    stop("No feasible transmission configurations.");


  // ------------------------------------------------------------
  // Sample lineage count proportional to total length
  // ------------------------------------------------------------

  NumericVector prob_lin(poss_lin.size());
  double total = 0.0;

  for (std::size_t i = 0; i < poss_lin.size(); ++i) {
    prob_lin[i] = len_lin[poss_lin[i]];
    total += prob_lin[i];
  }

  prob_lin = prob_lin / total;

  int idx = sample_index(prob_lin);   // returns 1-based index
  int lb = poss_lin[idx - 1];     // convert back to 1-based lb

  // ------------------------------------------------------------
  // Sample configuration within selected lineage count
  // ------------------------------------------------------------
  List configs = tr_lin[lb];
  const int K = configs.size();

  NumericVector len_conf(K);
  for (int k = 0; k < K; ++k)
    len_conf[k] = as<double>(as<List>(configs[k])["length"]);

  NumericVector prob_conf = len_conf / sum(len_conf);
  int config_idx = sample_index(prob_conf);

  List chosen = configs[config_idx - 1];

  // ------------------------------------------------------------
  // Sample transmission time uniformly over interval
  // ------------------------------------------------------------
  double start = as<double>(chosen["start"]);
  double end   = as<double>(chosen["end"]);
  double sam_time = R::runif(start, end);

  // ------------------------------------------------------------
  // Proposal density
  // ------------------------------------------------------------
  double prop_density =
    std::log(prob_lin[idx - 1]) +
    std::log(prob_conf[config_idx - 1]) -
    std::log(end - start);

  return List::create(
    _["lb"]           = lb,
    _["config_idx"]   = config_idx,
    _["branches"]     = chosen["branches"],
    _["bin"]          = chosen["bin"],
    _["sam_time"]     = sam_time,
    _["prop_density"] = prop_density
  );
}

