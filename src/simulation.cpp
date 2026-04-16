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
List sim_svp_norm_cpp(arma::vec beta, int p, int N, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);

  arma::mat veta = rnorm(n);
  arma::mat epst = rnorm(n);

  arma::mat v0(p, 1, arma::fill::zeros);
  h.rows(0, p - 1) = v0.rows(0, p - 1);
  u.row(0) = epst.row(0) * exp(h.row(0) / 2) * sy;
  for (int i = p; i < n; i++) {
    h.row(i) = v0(0) + trans(phi) * flipud(h.rows(i - p, i - 1)) + veta.row(i) * sv;
  }
  for (int i = 1; i < n; i++) {
    u.row(i) = epst.row(i) * exp(h.row(i) / 2) * sy;
  }
  u    = u.rows(burnin, n - 1);
  h    = h.rows(burnin, n - 1);
  epst = epst.rows(burnin, n - 1);
  veta = veta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["z"] = epst;
  out["v"] = veta;
  return(out);
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

  arma::mat v0(p, 1, arma::fill::zeros);
  h.rows(0, p - 1) = v0.rows(0, p - 1);
  u.row(0) = zeta.row(0) * exp(h.row(0) / 2) * sy;
  for (int i = p; i < n; i++) {
    h.row(i) = v0(0) + trans(phi) * flipud(h.rows(i - p, i - 1)) + veta.row(i) * sv;
  }
  for (int i = 1; i < n; i++) {
    u.row(i) = zeta.row(i) * exp(h.row(i) / 2) * sy;
  }
  u    = u.rows(burnin, n - 1);
  h    = h.rows(burnin, n - 1);
  zeta = zeta.rows(burnin, n - 1);
  veta = veta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["z"] = zeta;    // For Gaussian leverage, latent zeta == effective z (both N(0,1))
  out["v"] = veta;
  return(out);
}


// [[Rcpp::export]]
List sim_sv_t_cpp(arma::vec beta, int p, int N, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double nu = beta(p + 2);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);
  arma::mat eta = arma::randn(n, 1);
  arma::mat eps(n, 1);

  for (int i = 0; i < n; ++i) {
    eps(i, 0) = R::rt(nu);
  }

  for (int i = 0; i < p; ++i) {
    h(i, 0) = 0.0;
    u(i, 0) = eps(i, 0) * std::exp(h(i, 0) / 2.0) * sy;
  }
  for (int i = p; i < n; ++i) {
    h(i, 0) = arma::dot(phi, arma::flipud(h.col(0).subvec(i - p, i - 1))) + sv * eta(i, 0);
    u(i, 0) = eps(i, 0) * std::exp(h(i, 0) / 2.0) * sy;
  }

  u   = u.rows(burnin, n - 1);
  h   = h.rows(burnin, n - 1);
  eps = eps.rows(burnin, n - 1);
  eta = eta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["z"] = eps;   // t(nu)-distributed return innovation
  out["v"] = eta;   // N(0,1) volatility innovation
  return(out);
}


// [[Rcpp::export]]
List sim_sv_ged_cpp(arma::vec beta, int p, int N, int burnin) {
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double ged_nu = beta(p + 2);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);
  arma::mat eta(n, 1, arma::fill::randn);
  arma::mat eps = rged_cpp(n, 0, 1, ged_nu);
  eps.reshape(n, 1);

  for (int i = 0; i < p; ++i) {
    h(i, 0) = 0.0;
    u(i, 0) = eps(i, 0) * std::exp(h(i, 0) / 2.0) * sy;
  }
  for (int i = p; i < n; ++i) {
    h(i, 0) = arma::dot(phi, arma::flipud(h.col(0).subvec(i - p, i - 1))) + sv * eta(i, 0);
    u(i, 0) = eps(i, 0) * std::exp(h(i, 0) / 2.0) * sy;
  }

  u   = u.rows(burnin, n - 1);
  h   = h.rows(burnin, n - 1);
  eps = eps.rows(burnin, n - 1);
  eta = eta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["z"] = eps;   // unit-variance GED(nu) return innovation
  out["v"] = eta;   // N(0,1) volatility innovation
  return(out);
}


// Forward declaration from utils_cpp.cpp
double qged_std_cpp(double p, double nu);

// [[Rcpp::export]]
List sim_svp_leverage_t_cpp(arma::vec beta, int p, int N, int burnin) {
  // SVL(p) with Student-t errors via scale mixture
  // beta = (phi_1, ..., phi_p, sigma_y, sigma_v, nu, rho)
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double nu = beta(p + 2);
  double rho0 = beta(p + 3);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);
  arma::mat z_eff(n, 1, arma::fill::zeros);

  // Gaussian innovations for leverage structure
  arma::mat veta = rnorm(n);
  arma::mat eps = rnorm(n);
  arma::mat zeta = rnorm(n);  // latent Gaussian (drives leverage)

  // Correlation: v_t = rho * z_{t-1} + sqrt(1-rho^2) * eps_t
  veta.rows(1, n - 1) = rho0 * zeta.rows(0, n - 2) +
    std::sqrt(1.0 - std::pow(rho0, 2)) * eps.rows(1, n - 1);

  // Chi-squared mixing: nu*lambda ~ chi^2(nu), lambda = chi^2/nu
  arma::vec lambda_inv_sqrt(n);
  for (int i = 0; i < n; ++i) {
    double chi2_val = R::rchisq(nu);
    lambda_inv_sqrt(i) = std::sqrt(nu / chi2_val);  // lambda^{-1/2} = sqrt(nu/chi2)
  }

  // Log-volatility AR(p) recursion
  arma::mat v0(p, 1, arma::fill::zeros);
  h.rows(0, p - 1) = v0.rows(0, p - 1);
  z_eff(0, 0) = zeta(0, 0) * lambda_inv_sqrt(0);
  u(0, 0) = z_eff(0, 0) * std::exp(h(0, 0) / 2.0) * sy;
  for (int i = p; i < n; i++) {
    h(i, 0) = arma::dot(phi, arma::flipud(h.col(0).subvec(i - p, i - 1))) + veta(i, 0) * sv;
  }
  for (int i = 1; i < n; i++) {
    z_eff(i, 0) = zeta(i, 0) * lambda_inv_sqrt(i);
    u(i, 0) = z_eff(i, 0) * std::exp(h(i, 0) / 2.0) * sy;
  }

  u     = u.rows(burnin, n - 1);
  h     = h.rows(burnin, n - 1);
  z_eff = z_eff.rows(burnin, n - 1);
  veta  = veta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["z"] = z_eff;   // effective t(nu) return innovation: z = zeta * lambda^{-1/2}
  out["v"] = veta;    // N(0,1) volatility innovation (includes leverage correlation)
  return(out);
}


// [[Rcpp::export]]
List sim_svp_leverage_ged_cpp(arma::vec beta, int p, int N, int burnin) {
  // SVL(p) with GED errors via Gaussian copula
  // beta = (phi_1, ..., phi_p, sigma_y, sigma_v, nu_ged, rho)
  int n = burnin + N;
  arma::vec phi = beta.subvec(0, p - 1);
  double sy = beta(p);
  double sv = beta(p + 1);
  double ged_nu = beta(p + 2);
  double rho0 = beta(p + 3);

  arma::mat h(n, 1, arma::fill::zeros);
  arma::mat u(n, 1, arma::fill::zeros);
  arma::mat z_eff(n, 1, arma::fill::zeros);

  // Gaussian innovations for leverage structure
  arma::mat veta = rnorm(n);
  arma::mat eps = rnorm(n);
  arma::mat zeta = rnorm(n);  // latent Gaussian

  // Correlation: v_t = rho * z_{t-1} + sqrt(1-rho^2) * eps_t
  veta.rows(1, n - 1) = rho0 * zeta.rows(0, n - 2) +
    std::sqrt(1.0 - std::pow(rho0, 2)) * eps.rows(1, n - 1);

  // Log-volatility AR(p) recursion
  arma::mat v0(p, 1, arma::fill::zeros);
  h.rows(0, p - 1) = v0.rows(0, p - 1);
  // u_t = sigma_y * exp(h_t/2) * F_GED^{-1}(Phi(z_t))
  double p_t = R::pnorm(zeta(0, 0), 0.0, 1.0, 1, 0);
  z_eff(0, 0) = qged_std_cpp(p_t, ged_nu);
  u(0, 0) = z_eff(0, 0) * std::exp(h(0, 0) / 2.0) * sy;
  for (int i = p; i < n; i++) {
    h(i, 0) = arma::dot(phi, arma::flipud(h.col(0).subvec(i - p, i - 1))) + veta(i, 0) * sv;
  }
  for (int i = 1; i < n; i++) {
    p_t = R::pnorm(zeta(i, 0), 0.0, 1.0, 1, 0);
    z_eff(i, 0) = qged_std_cpp(p_t, ged_nu);
    u(i, 0) = z_eff(i, 0) * std::exp(h(i, 0) / 2.0) * sy;
  }

  u     = u.rows(burnin, n - 1);
  h     = h.rows(burnin, n - 1);
  z_eff = z_eff.rows(burnin, n - 1);
  veta  = veta.rows(burnin, n - 1);

  List out;
  out["y"] = u;
  out["h"] = h;
  out["z"] = z_eff;   // effective GED(nu) return innovation: z = F_GED^{-1}(Phi(zeta))
  out["v"] = veta;    // N(0,1) volatility innovation (includes leverage correlation)
  return(out);
}
