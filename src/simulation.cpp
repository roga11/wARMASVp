// ============================================================================ //
// wARMASVp: Simulation Functions for SV(p) Models
// Authors: Nazmul Ahsan, Jean-Marie Dufour, Gabriel Rodriguez Rondon
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declaration from utils_cpp.cpp
arma::vec rged_cpp(int n, double mean, double sd, double nu);

// [[Rcpp::export]]
arma::mat sim_svp_norm_cpp(arma::vec beta, int p, int N, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);

  arma::mat veta = rnorm(n + 1);
  arma::mat epst = rnorm(n + 1);

  arma::mat v(p, 1, arma::fill::zeros);
  h.rows(0, p - 1) = v.rows(0, p - 1);
  u.row(0) = epst.row(0) * exp(h.row(0) / 2) * sy;
  for (int i = p; i < n; i++) {
    h.row(i) = v(0) + trans(phi) * flipud(h.rows(i - p, i - 1)) + veta.row(i) * sv;
  }
  for (int i = 1; i < n; i++) {
    u.row(i) = epst.row(i - 1) * exp(h.row(i) / 2) * sy;
  }
  u = u.rows(burnin, n - 1);
  return(u);
}


// [[Rcpp::export]]
List sim_svp_leverage_norm_cpp(arma::vec beta, int p, int N, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double rho0 = beta(p + 2);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);

  arma::mat veta = rnorm(n);
  arma::mat eps = rnorm(n);
  arma::mat zeta = rnorm(n);

  veta.rows(1, n - 1) = rho0 * zeta.rows(0, n - 2) + sqrt(1 - pow(rho0, 2)) * eps.rows(1, n - 1);

  arma::mat v(p, 1, arma::fill::zeros);
  h.rows(0, p - 1) = v.rows(0, p - 1);
  u.row(0) = zeta.row(0) * exp(h.row(0) / 2) * sy;
  for (int i = p; i < n; i++) {
    h.row(i) = v(0) + trans(phi) * flipud(h.rows(i - p, i - 1)) + veta.row(i) * sv;
  }
  for (int i = 1; i < n; i++) {
    u.row(i) = zeta.row(i) * exp(h.row(i) / 2) * sy;
  }
  u = u.rows(burnin, n - 1);
  h = h.rows(burnin, n - 1);
  zeta = zeta.rows(burnin, n - 1);
  veta = veta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["zeta"] = zeta;
  out["veta"] = veta;
  return(out);
}


// [[Rcpp::export]]
arma::mat sim_sv_t_cpp(arma::vec beta, int p, int N, int K, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double nu = beta(p + 2);

  arma::mat h(n, K, arma::fill::zeros);
  arma::mat u(n, K, arma::fill::zeros);
  arma::mat eta = arma::randn(n, K);
  arma::mat eps(n, K);

  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < n; ++i) {
      eps(i, j) = R::rt(nu);
    }
  }

  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < p; ++i) {
      h(i, j) = 0.0;
      u(i, j) = eps(i, j) * std::exp(h(i, j) / 2.0) * sy;
    }
    for (int i = p; i < n; ++i) {
      h(i, j) = arma::dot(phi, arma::flipud(h.col(j).subvec(i - p, i - 1))) + sv * eta(i, j);
      u(i, j) = eps(i, j) * std::exp(h(i, j) / 2.0) * sy;
    }
  }
  return u.rows(burnin, n - 1);
}


// [[Rcpp::export]]
arma::mat sim_sv_ged_cpp(arma::vec beta, int p, int N, int K, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double ged_nu = beta(p + 2);

  arma::mat h(n, K, arma::fill::zeros);
  arma::mat u(n, K, arma::fill::zeros);
  arma::mat eta(n, K, arma::fill::randn);
  arma::mat eps = rged_cpp(n * K, 0, 1, ged_nu);
  eps.reshape(n, K);

  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < p; ++i) {
      h(i, j) = 0.0;
      u(i, j) = eps(i, j) * std::exp(h(i, j) / 2.0) * sy;
    }
    for (int i = p; i < n; ++i) {
      h(i, j) = arma::dot(phi, arma::flipud(h.col(j).subvec(i - p, i - 1))) + sv * eta(i, j);
      u(i, j) = eps(i, j) * std::exp(h(i, j) / 2.0) * sy;
    }
  }
  u = u.rows(burnin, n - 1);
  return(u);
}
