// ============================================================================ //
// wARMASVp: Simulation from Pre-Drawn Innovations (Fixed-Error MMC)
//
// These functions reconstruct time series y from pre-drawn innovation vectors
// and model parameters. Used for the MMC procedure per Dufour (2006, eq 4.22):
// fix disturbance vectors, vary only nuisance parameters.
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declaration from utils_cpp.cpp
double qged_std_cpp(double p, double nu);

// --------------------------------------------------------------------------- //
// Gaussian, no leverage
// Innovations: eta_vec (state, N(0,1)), eps_vec (measurement, N(0,1))
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::vec sim_from_innov_gaussian_cpp(const arma::vec& phi,
                                      double sigma_y, double sigma_v,
                                      const arma::vec& eta_vec,
                                      const arma::vec& eps_vec,
                                      int p, int T_out, int burnin) {
  int n = burnin + T_out;
  arma::vec h(n, arma::fill::zeros);
  arma::vec y(n, arma::fill::zeros);

  y(0) = eps_vec(0) * std::exp(h(0) / 2.0) * sigma_y;
  for (int i = p; i < n; i++) {
    h(i) = arma::dot(phi, arma::flipud(h.subvec(i - p, i - 1))) + sigma_v * eta_vec(i);
  }
  for (int i = 1; i < n; i++) {
    y(i) = eps_vec(i) * std::exp(h(i) / 2.0) * sigma_y;
  }
  return y.subvec(burnin, n - 1);
}

// --------------------------------------------------------------------------- //
// Gaussian, leverage
// Innovations: zeta_vec (N(0,1)), aux_vec (N(0,1))
// veta is derived: veta_t = rho*zeta_{t-1} + sqrt(1-rho^2)*aux_t
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::vec sim_from_innov_gaussian_lev_cpp(const arma::vec& phi,
                                          double sigma_y, double sigma_v,
                                          double rho,
                                          const arma::vec& zeta_vec,
                                          const arma::vec& aux_vec,
                                          int p, int T_out, int burnin) {
  int n = burnin + T_out;
  arma::vec h(n, arma::fill::zeros);
  arma::vec y(n, arma::fill::zeros);
  arma::vec veta(n, arma::fill::zeros);

  double sqrt_1mrho2 = std::sqrt(1.0 - rho * rho);
  for (int i = 1; i < n; i++) {
    veta(i) = rho * zeta_vec(i - 1) + sqrt_1mrho2 * aux_vec(i);
  }

  y(0) = zeta_vec(0) * std::exp(h(0) / 2.0) * sigma_y;
  for (int i = p; i < n; i++) {
    h(i) = arma::dot(phi, arma::flipud(h.subvec(i - p, i - 1))) + sigma_v * veta(i);
  }
  for (int i = 1; i < n; i++) {
    y(i) = zeta_vec(i) * std::exp(h(i) / 2.0) * sigma_y;
  }
  return y.subvec(burnin, n - 1);
}

// --------------------------------------------------------------------------- //
// Student-t, no leverage
// Innovations: eta_vec (state, N(0,1)), eps_vec (measurement, pre-drawn t(nu))
// nu is fixed (not a nuisance parameter for this test)
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::vec sim_from_innov_t_cpp(const arma::vec& phi,
                               double sigma_y, double sigma_v,
                               const arma::vec& eta_vec,
                               const arma::vec& eps_vec,
                               int p, int T_out, int burnin) {
  int n = burnin + T_out;
  arma::vec h(n, arma::fill::zeros);
  arma::vec y(n, arma::fill::zeros);

  for (int i = 0; i < p; i++) {
    y(i) = eps_vec(i) * std::exp(h(i) / 2.0) * sigma_y;
  }
  for (int i = p; i < n; i++) {
    h(i) = arma::dot(phi, arma::flipud(h.subvec(i - p, i - 1))) + sigma_v * eta_vec(i);
    y(i) = eps_vec(i) * std::exp(h(i) / 2.0) * sigma_y;
  }
  return y.subvec(burnin, n - 1);
}

// --------------------------------------------------------------------------- //
// GED, no leverage
// Innovations: eta_vec (state, N(0,1)), eps_vec (measurement, pre-drawn GED(nu))
// nu is fixed (not a nuisance parameter for this test)
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::vec sim_from_innov_ged_cpp(const arma::vec& phi,
                                 double sigma_y, double sigma_v,
                                 const arma::vec& eta_vec,
                                 const arma::vec& eps_vec,
                                 int p, int T_out, int burnin) {
  // Same structure as Student-t — eps already drawn from GED
  return sim_from_innov_t_cpp(phi, sigma_y, sigma_v, eta_vec, eps_vec,
                              p, T_out, burnin);
}

// --------------------------------------------------------------------------- //
// Student-t, leverage
// Innovations: zeta_vec (N(0,1)), aux_vec (N(0,1)), U_chi2_vec (Uniform(0,1))
// nu VARIES (nuisance parameter) — chi2 reconstructed via PIT: qchisq(U, nu)
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::vec sim_from_innov_t_lev_cpp(const arma::vec& phi,
                                   double sigma_y, double sigma_v,
                                   double nu, double rho,
                                   const arma::vec& zeta_vec,
                                   const arma::vec& aux_vec,
                                   const arma::vec& U_chi2_vec,
                                   int p, int T_out, int burnin) {
  int n = burnin + T_out;
  arma::vec h(n, arma::fill::zeros);
  arma::vec y(n, arma::fill::zeros);
  arma::vec veta(n, arma::fill::zeros);

  double sqrt_1mrho2 = std::sqrt(1.0 - rho * rho);
  for (int i = 1; i < n; i++) {
    veta(i) = rho * zeta_vec(i - 1) + sqrt_1mrho2 * aux_vec(i);
  }

  // Chi-squared mixing via PIT: chi2 = qchisq(U, nu)
  // Guard against U=0 (machine zero) which would give chi2=0 and sqrt(nu/0)=Inf
  arma::vec lambda_inv_sqrt(n);
  for (int i = 0; i < n; i++) {
    double chi2_val = R::qchisq(U_chi2_vec(i), nu, 1, 0);
    chi2_val = std::max(chi2_val, 1e-300);
    lambda_inv_sqrt(i) = std::sqrt(nu / chi2_val);
  }

  y(0) = zeta_vec(0) * lambda_inv_sqrt(0) * std::exp(h(0) / 2.0) * sigma_y;
  for (int i = p; i < n; i++) {
    h(i) = arma::dot(phi, arma::flipud(h.subvec(i - p, i - 1))) + sigma_v * veta(i);
  }
  for (int i = 1; i < n; i++) {
    y(i) = zeta_vec(i) * lambda_inv_sqrt(i) * std::exp(h(i) / 2.0) * sigma_y;
  }
  return y.subvec(burnin, n - 1);
}

// --------------------------------------------------------------------------- //
// GED, leverage
// Innovations: zeta_vec (N(0,1)), aux_vec (N(0,1))
// nu VARIES — GED measurement via Gaussian copula: u_t = qged(pnorm(zeta_t), nu)
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::vec sim_from_innov_ged_lev_cpp(const arma::vec& phi,
                                     double sigma_y, double sigma_v,
                                     double nu, double rho,
                                     const arma::vec& zeta_vec,
                                     const arma::vec& aux_vec,
                                     int p, int T_out, int burnin) {
  int n = burnin + T_out;
  arma::vec h(n, arma::fill::zeros);
  arma::vec y(n, arma::fill::zeros);
  arma::vec veta(n, arma::fill::zeros);

  double sqrt_1mrho2 = std::sqrt(1.0 - rho * rho);
  for (int i = 1; i < n; i++) {
    veta(i) = rho * zeta_vec(i - 1) + sqrt_1mrho2 * aux_vec(i);
  }

  // GED via Gaussian copula: u_t = F_GED^{-1}(Phi(zeta_t))
  double p_t = R::pnorm(zeta_vec(0), 0.0, 1.0, 1, 0);
  y(0) = qged_std_cpp(p_t, nu) * std::exp(h(0) / 2.0) * sigma_y;
  for (int i = p; i < n; i++) {
    h(i) = arma::dot(phi, arma::flipud(h.subvec(i - p, i - 1))) + sigma_v * veta(i);
  }
  for (int i = 1; i < n; i++) {
    p_t = R::pnorm(zeta_vec(i), 0.0, 1.0, 1, 0);
    y(i) = qged_std_cpp(p_t, nu) * std::exp(h(i) / 2.0) * sigma_y;
  }
  return y.subvec(burnin, n - 1);
}
