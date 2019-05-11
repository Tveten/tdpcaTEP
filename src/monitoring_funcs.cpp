#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix update_cumsum(NumericMatrix & u, NumericVector next_x, int t) {
  u(_, t) = u(_, t - 1) + next_x;
  return u;
}

// [[Rcpp::export]]
NumericMatrix update_cumsum2(NumericMatrix & v,
                             NumericMatrix & u,
                             NumericVector next_x, int m, int t) {
  NumericVector dev2 = pow(next_x - u(_, t - 1) / (m + t - 1), 2);
  v(_, t) = v(_, t - 1) + (m + t - 1.0) / (m + t) * dev2;
  return v;
}

// [[Rcpp::export]]
List update_sumsC(List & sums, NumericVector next_x, int m, int t) {
  NumericMatrix u = as<NumericMatrix>(sums["u"]);
  NumericMatrix v = as<NumericMatrix>(sums["v"]);

  u = update_cumsum(u, next_x, t);
  v = update_cumsum2(v, u, next_x, m, t);

  return sums;
}

// [[Rcpp::export]]
NumericVector log_lik_ratio(int m, int t, int k,
                            NumericVector& s_full,
                            NumericVector& s_pre,
                            NumericVector& s_post) {
  if (is_true(any(s_post <= 0))) {
    NumericVector log_lik(s_full.size());
    return log_lik;
  } else {
    NumericVector log_lik = (m + t) * log(s_full) -
                            (m + k) * log(s_pre) -
                            (t - k) * log(s_post);
    return log_lik;
  }
}

// [[Rcpp::export]]
NumericVector log_liks_matC(List & sums, int m, int t, int w, double p0) {
  NumericMatrix u = as<NumericMatrix>(sums["u"]);
  NumericMatrix v = as<NumericMatrix>(sums["v"]);

  NumericVector u_full = u(_, t);
  NumericVector v_full = v(_, t);
  NumericVector s_full = v_full / (m + t);
  int d = u_full.size();

  int k_start = std::max(0, t - (w + 1));
  int k_end   = t - 2;
  NumericMatrix log_liks(d, k_end - k_start + 1);
  for(int k = k_start; k <= k_end; ++k) {
    NumericVector u_pre   = u(_, k);
    NumericVector v_pre   = v(_, k);
    NumericVector s_pre   = v_pre / (m + k);
    NumericVector avg_pre = u_pre / (m + k);

    NumericVector avg_post = (u_full - u_pre) / (t - k);

    NumericVector v_post   = v_full - v_pre - ((double)(m + k) * (t - k)) / (m + t) *
                               pow((avg_pre - avg_post), 2);
    NumericVector s_post   = 1.0 / (t - k) * v_post;
    int i = k - k_start;
    log_liks(_, i) = log_lik_ratio(m, t, k, s_full, s_pre, s_post);
  }
  return log_liks;
}
