// ============================================================================ //
// wARMASVp: GMM Moment Functions for Hypothesis Testing
// Authors: Nazmul Ahsan, Jean-Marie Dufour, Gabriel Rodriguez Rondon
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declarations from utils_cpp.cpp
double kendall_corr(arma::vec x, arma::vec y);
double qged_std_cpp(double p, double nu);

// --------------------------------------------------------------------------- //
// Internal helpers
// --------------------------------------------------------------------------- //

// Batch autocovariance: compute gamma(0), gamma(1), ..., gamma(max_lag)
// Each lag k uses normalization 1/(T-k), matching acov_g() in utils_cpp.cpp
static arma::vec batch_acov(const arma::vec& ys, int max_lag) {
  int T = ys.n_elem;
  arma::vec gamk(max_lag + 1);
  for (int k = 0; k <= max_lag; k++) {
    double s = 0.0;
    for (int t = k; t < T; t++) {
      s += ys(t) * ys(t - k);
    }
    gamk(k) = s / (T - k);
  }
  return gamk;
}

// Per-observation lagged products: gamk_t(t, k) = ys(t) * ys(t-k) for k=1..max_lag
// gamk_t has T rows, max_lag columns. First k rows of column k are 0.
static arma::mat build_gamk_t(const arma::vec& ys, int max_lag) {
  int T = ys.n_elem;
  arma::mat gamk_t(T, max_lag, arma::fill::zeros);
  for (int k = 0; k < max_lag; k++) {
    int lag = k + 1;
    for (int t = lag; t < T; t++) {
      gamk_t(t, k) = ys(t) * ys(t - lag);
    }
  }
  return gamk_t;
}

// HAC Bartlett kernel covariance estimation
// g_t: T_eff x n_mom matrix of per-observation scores
// Tsize: original sample size (used for bandwidth)
// Returns W (n_mom x n_mom) or Gam0 if Bartlett=false
// pinv_fn: R function for pseudoinverse
static arma::mat hac_estimate(const arma::mat& g_t, int Tsize, bool Bartlett,
                               Rcpp::Function pinv_fn) {
  int T_eff = g_t.n_rows;
  arma::mat Gam0 = (1.0 / T_eff) * (g_t.t() * g_t);
  if (Bartlett) {
    int KT = (int)std::floor(std::pow((double)Tsize, 1.0 / 3.0));
    arma::mat W = Gam0;
    for (int k = 1; k <= KT; k++) {
      double weight = 1.0 - ((double)k / (KT + 1.0));
      arma::mat Gamk = g_t.rows(k, T_eff - 1).t() * g_t.rows(0, T_eff - 1 - k) / (T_eff - 1.0);
      W += weight * (Gamk + Gamk.t());
    }
    return Rcpp::as<arma::mat>(pinv_fn(Rcpp::wrap(W)));
  } else {
    return Rcpp::as<arma::mat>(pinv_fn(Rcpp::wrap(Gam0)));
  }
}

// Compute EH (leverage cross-moment) — Pearson or Kendall
// y: T x 1 column vector (original returns)
// Returns scalar EH
static double compute_EH(const arma::vec& y, const std::string& rho_type) {
  int T = y.n_elem;
  arma::vec yabs = arma::abs(y);
  double muu = arma::mean(y.subvec(0, T - 2));
  double mua = arma::mean(yabs.subvec(1, T - 1));
  if (rho_type == "kendall") {
    arma::vec x1 = yabs.subvec(1, T - 1) - mua;
    arma::vec x2 = y.subvec(0, T - 2) - muu;
    double tau = kendall_corr(x1, x2);
    // var using N-1 denominator (stats::var)
    double var1 = arma::as_scalar(arma::var(yabs.subvec(1, T - 1)));
    double var2 = arma::as_scalar(arma::var(y.subvec(0, T - 2)));
    return tau * std::sqrt(var1 * var2);
  } else {
    arma::vec x1 = yabs.subvec(1, T - 1) - mua;
    arma::vec x2 = y.subvec(0, T - 2) - muu;
    return arma::dot(x1, x2) / (T - 2.0);
  }
}

// Compute gamtmp = sum_{j=2}^{p} phi[j] * (gamk[j-1] + gamk[j])
// gamk is 0-indexed: gamk(0)=gamma(1), gamk(1)=gamma(2), etc.
static double compute_gamtmp(const arma::vec& phi, const arma::vec& gamk, int p) {
  double gamtmp = 0.0;
  for (int j = 1; j < p; j++) { // j=1 means phi[2] in 1-based
    gamtmp += phi(j) * (gamk(j - 1) + gamk(j)); // gamk(j-1) = gamma(j), gamk(j) = gamma(j+1)
  }
  return gamtmp;
}

// Compute per-obs gamtmp_t = sum_{j=2}^{p} phi[j] * (gamk_t[,j-1] + gamk_t[,j])
static arma::vec compute_gamtmp_t(const arma::vec& phi, const arma::mat& gamk_t, int p) {
  int T = gamk_t.n_rows;
  arma::vec gamtmp_t(T, arma::fill::zeros);
  for (int j = 1; j < p; j++) {
    gamtmp_t += phi(j) * (gamk_t.col(j - 1) + gamk_t.col(j));
  }
  return gamtmp_t;
}

// Build YW moment conditions: mk[j] = gamk[p+j] - sum(phi * gamk[p+j-(1:p)])
// gamk_all is 0-indexed with gamk_all(0)=gamma(0), gamk_all(1)=gamma(1), etc.
static arma::vec build_mk(const arma::vec& phi, const arma::vec& gamk_all, int p) {
  arma::vec mk(p);
  for (int j = 0; j < p; j++) {
    int lag_idx = p + j + 1; // gamma(p+j+1) in 1-based = gamk_all(p+j+1)
    double s = 0.0;
    for (int i = 0; i < p; i++) {
      s += phi(i) * gamk_all(lag_idx - i - 1);
    }
    mk(j) = gamk_all(lag_idx) - s;
  }
  return mk;
}

// Build per-obs YW score columns
static arma::mat build_mk_t(const arma::vec& phi, const arma::mat& gamk_t, int p) {
  int T = gamk_t.n_rows;
  arma::mat mk_t(T, p);
  for (int j = 0; j < p; j++) {
    int lag_idx = p + j; // 0-based column in gamk_t for gamma(p+j+1)
    arma::vec col = gamk_t.col(lag_idx);
    for (int i = 0; i < p; i++) {
      col -= phi(i) * gamk_t.col(lag_idx - i - 1);
    }
    mk_t.col(j) = col;
  }
  return mk_t;
}


// =========================================================================== //
// Fixed-Amat functions (pure C++, no .pinv needed)
// =========================================================================== //

// [[Rcpp::export]]
double LRT_moment_lev_svp_cpp(arma::vec y, Rcpp::List mdl_out,
                                arma::mat Amat, std::string rho_type,
                                double del = 1e-10) {
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double rho = mdl_out["rho"];
  double gammatilde = mdl_out["gammatilde"];

  // Batch autocovariances: gamma(0) through gamma(2p)
  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);
  // gamk(1..2p) for compatibility: gamk[xp] = gamk_all(xp) for xp=1..2p
  double gamtmp = compute_gamtmp(phi, gamk_all.subvec(1, 2 * p), p);

  // EH
  double EH = compute_EH(y, rho_type);

  // Moment vector
  double m1 = mu + 1.2704 - std::log(sigy * sigy);
  double m2 = gam0 + gamk_all(1) - (M_PI * M_PI / 2.0) -
    (1.0 / (1.0 - phi(0))) * (gamtmp + sigv * sigv);
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m4 = rho - (EH * std::sqrt(2.0 * M_PI)) /
    (sigv * sigy * sigy) * std::exp(-0.25 * gammatilde);

  arma::vec g(p + 3);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m4;

  return arma::as_scalar(g.t() * Amat * g);
}


// [[Rcpp::export]]
double LRT_moment_lev_t_cpp(arma::vec y, Rcpp::List mdl_out,
                             arma::mat Amat, std::string rho_type,
                             double del = 1e-10) {
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double rho = mdl_out["rho"];
  double gammatilde = mdl_out["gammatilde"];
  double nu = mdl_out["v"];

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);

  double EH = compute_EH(y, rho_type);

  // Theoretical lag-1 ACF via ARMAacf
  Rcpp::Function ARMAacf("ARMAacf", Rcpp::Environment::namespace_env("stats"));
  arma::vec acf_vals = Rcpp::as<arma::vec>(ARMAacf(Rcpp::Named("ar") = phi,
                                                     Rcpp::Named("lag.max") = 1));
  double rho_w1 = acf_vals(1);

  // Distribution-specific constants
  double mu_log_sq_t = R::psigamma(0.5, 0) - R::psigamma(nu / 2.0, 0) + std::log(nu);
  double var_log_sq_t = R::psigamma(0.5, 1) + R::psigamma(nu / 2.0, 1);
  // C_t(nu) = (nu/2) * [Gamma((nu-1)/2) / Gamma(nu/2)]^2
  double Ct = std::exp(std::log(nu / 2.0) +
    2.0 * (R::lgammafn((nu - 1.0) / 2.0) - R::lgammafn(nu / 2.0)));

  double m1 = mu - mu_log_sq_t - std::log(sigy * sigy);
  double m2 = gam0 - arma::dot(phi, gamk_all.subvec(1, p)) - var_log_sq_t - sigv * sigv;
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m_prof = var_log_sq_t - gam0 + gamk_all(1) / rho_w1;
  double m_lev = rho - (std::sqrt(2.0 * M_PI) * EH) /
    (sigv * sigy * sigy * Ct) * std::exp(-0.25 * gammatilde);

  arma::vec g(p + 4);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m_prof;
  g(p + 3) = m_lev;

  return arma::as_scalar(g.t() * Amat * g);
}


// [[Rcpp::export]]
double LRT_moment_lev_ged_cpp(arma::vec y, Rcpp::List mdl_out,
                               arma::mat Amat, std::string rho_type,
                               double del,
                               arma::vec gh_nodes,
                               arma::vec gh_weights) {
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double rho = mdl_out["rho"];
  double gammatilde = mdl_out["gammatilde"];
  double nu = mdl_out["v"];

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);

  double EH = compute_EH(y, rho_type);

  Rcpp::Function ARMAacf("ARMAacf", Rcpp::Environment::namespace_env("stats"));
  arma::vec acf_vals = Rcpp::as<arma::vec>(ARMAacf(Rcpp::Named("ar") = phi,
                                                     Rcpp::Named("lag.max") = 1));
  double rho_w1 = acf_vals(1);

  // GED constants
  double mu_log_sq_ged = (2.0 / nu) * R::psigamma(1.0 / nu, 0) +
    R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu);
  double var_log_sq_ged = std::pow(2.0 / nu, 2) * R::psigamma(1.0 / nu, 1);
  // E[|u|] for standardized GED
  double E_abs_u = std::exp(0.5 * (R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu)) +
    R::lgammafn(2.0 / nu) - R::lgammafn(1.0 / nu));

  // GED leverage moment via Gauss-Hermite quadrature
  double c_shift = rho * sigv / 2.0;
  double E_g_shifted = 0.0;
  for (int i = 0; i < (int)gh_nodes.n_elem; i++) {
    double p_val = R::pnorm5(gh_nodes(i) + c_shift, 0.0, 1.0, 1, 0);
    E_g_shifted += gh_weights(i) * qged_std_cpp(p_val, nu);
  }

  double m1 = mu - mu_log_sq_ged - std::log(sigy * sigy);
  double m2 = gam0 - arma::dot(phi, gamk_all.subvec(1, p)) - var_log_sq_ged - sigv * sigv;
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m_prof = var_log_sq_ged - gam0 + gamk_all(1) / rho_w1;
  double m_lev = sigy * sigy * E_abs_u * E_g_shifted *
    std::exp(gammatilde / 4.0) - EH;

  arma::vec g(p + 4);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m_prof;
  g(p + 3) = m_lev;

  return arma::as_scalar(g.t() * Amat * g);
}


// [[Rcpp::export]]
double LRT_moment_t_cpp(arma::vec y, Rcpp::List mdl_out,
                         arma::mat Amat, bool WAmat = false,
                         double del = 1e-10, bool Bartlett = true,
                         Rcpp::Nullable<Rcpp::Function> pinv_fn = R_NilValue) {
  int Tsize = y.n_elem;
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double nu = mdl_out["v"];
  int n_mom = p + 3;

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);

  Rcpp::Function ARMAacf("ARMAacf", Rcpp::Environment::namespace_env("stats"));
  arma::vec acf_vals = Rcpp::as<arma::vec>(ARMAacf(Rcpp::Named("ar") = phi,
                                                     Rcpp::Named("lag.max") = 1));
  double rho_w1 = acf_vals(1);

  double mu_log_sq_t = R::psigamma(0.5, 0) - R::psigamma(nu / 2.0, 0) + std::log(nu);
  double var_log_sq_t = R::psigamma(0.5, 1) + R::psigamma(nu / 2.0, 1);

  double m1 = mu - mu_log_sq_t - std::log(sigy * sigy);
  double m2 = gam0 - arma::dot(phi, gamk_all.subvec(1, p)) - var_log_sq_t - sigv * sigv;
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m_prof = var_log_sq_t - gam0 + gamk_all(1) / rho_w1;

  arma::vec g(n_mom);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m_prof;

  arma::mat Amat_use = Amat;

  if (WAmat) {
    // Build per-observation score matrix
    arma::mat gamk_t = build_gamk_t(ys, 2 * p);
    arma::vec gam0_t = ys % ys;

    arma::mat g_t(Tsize, n_mom, arma::fill::zeros);
    g_t.col(0) = ly2 - mu_log_sq_t - std::log(sigy * sigy);
    // g2: gam0_t - sum(phi_j * gamk_t[,j]) - var_log_sq_t - sigv^2
    arma::vec g2 = gam0_t;
    for (int j = 0; j < p; j++) {
      g2 -= phi(j) * gamk_t.col(j);
    }
    g_t.col(1) = g2 - var_log_sq_t - sigv * sigv;

    // YW per-obs
    arma::mat mk_t_mat = build_mk_t(phi, gamk_t, p);
    g_t.cols(2, p + 1) = mk_t_mat;

    // Profiling per-obs
    g_t.col(p + 2) = var_log_sq_t - gam0_t + gamk_t.col(0) / rho_w1;

    // Trim initial p rows
    g_t = g_t.rows(p, Tsize - 1);

    // HAC
    Rcpp::Function pinv_f(pinv_fn.get());
    Amat_use = hac_estimate(g_t, Tsize, Bartlett, pinv_f);
  }

  return arma::as_scalar(g.t() * Amat_use * g);
}


// [[Rcpp::export]]
double LRT_moment_ged_cpp(arma::vec y, Rcpp::List mdl_out,
                           arma::mat Amat, bool WAmat = false,
                           double del = 1e-10, bool Bartlett = true,
                           Rcpp::Nullable<Rcpp::Function> pinv_fn = R_NilValue) {
  int Tsize = y.n_elem;
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double nu = mdl_out["v"];
  int n_mom = p + 3;

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);

  Rcpp::Function ARMAacf("ARMAacf", Rcpp::Environment::namespace_env("stats"));
  arma::vec acf_vals = Rcpp::as<arma::vec>(ARMAacf(Rcpp::Named("ar") = phi,
                                                     Rcpp::Named("lag.max") = 1));
  double rho_w1 = acf_vals(1);

  double mu_log_sq_ged = (2.0 / nu) * R::psigamma(1.0 / nu, 0) +
    R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu);
  double var_log_sq_ged = std::pow(2.0 / nu, 2) * R::psigamma(1.0 / nu, 1);

  double m1 = mu - mu_log_sq_ged - std::log(sigy * sigy);
  double m2 = gam0 - arma::dot(phi, gamk_all.subvec(1, p)) - var_log_sq_ged - sigv * sigv;
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m_prof = var_log_sq_ged - gam0 + gamk_all(1) / rho_w1;

  arma::vec g(n_mom);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m_prof;

  arma::mat Amat_use = Amat;

  if (WAmat) {
    arma::mat gamk_t = build_gamk_t(ys, 2 * p);
    arma::vec gam0_t = ys % ys;

    arma::mat g_t(Tsize, n_mom, arma::fill::zeros);
    g_t.col(0) = ly2 - mu_log_sq_ged - std::log(sigy * sigy);
    arma::vec g2 = gam0_t;
    for (int j = 0; j < p; j++) {
      g2 -= phi(j) * gamk_t.col(j);
    }
    g_t.col(1) = g2 - var_log_sq_ged - sigv * sigv;

    arma::mat mk_t_mat = build_mk_t(phi, gamk_t, p);
    g_t.cols(2, p + 1) = mk_t_mat;

    g_t.col(p + 2) = var_log_sq_ged - gam0_t + gamk_t.col(0) / rho_w1;

    g_t = g_t.rows(p, Tsize - 1);

    Rcpp::Function pinv_f(pinv_fn.get());
    Amat_use = hac_estimate(g_t, Tsize, Bartlett, pinv_f);
  }

  return arma::as_scalar(g.t() * Amat_use * g);
}


// =========================================================================== //
// _Amat functions (HAC weighting, call R's .pinv via passed function)
// =========================================================================== //

// [[Rcpp::export]]
double LRT_moment_lev_svp_Amat_cpp(arma::vec y, Rcpp::List mdl_out,
                                     std::string rho_type,
                                     double del = 1e-10, bool Bartlett = true,
                                     Rcpp::Nullable<Rcpp::Function> pinv_fn = R_NilValue) {
  int Tsize = y.n_elem;
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double rho = mdl_out["rho"];
  double gammatilde = mdl_out["gammatilde"];
  int n_mom = p + 3;

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);
  double gamtmp = compute_gamtmp(phi, gamk_all.subvec(1, 2 * p), p);

  // Per-obs products
  arma::mat gamk_t = build_gamk_t(ys, 2 * p);
  arma::vec gam0_t = ys % ys;
  arma::vec gamtmp_t = compute_gamtmp_t(phi, gamk_t, p);

  // EH + per-obs leverage
  arma::vec yabs = arma::abs(y);
  double muu = arma::mean(y.subvec(0, Tsize - 2));
  double mua = arma::mean(yabs.subvec(1, Tsize - 1));
  double EH;
  arma::mat g_t(Tsize, n_mom, arma::fill::zeros);

  if (rho_type == "kendall") {
    arma::vec x1 = yabs.subvec(1, Tsize - 1) - mua;
    arma::vec x2 = y.subvec(0, Tsize - 2) - muu;
    double tau = kendall_corr(x1, x2);
    double var1 = arma::as_scalar(arma::var(yabs.subvec(1, Tsize - 1)));
    double var2 = arma::as_scalar(arma::var(y.subvec(0, Tsize - 2)));
    EH = tau * std::sqrt(var1 * var2);
    // Per-obs: kendall_corr * sqrt(|y_{t+1}|^2 * y_t^2)
    for (int t = 1; t < Tsize; t++) {
      g_t(t, n_mom - 1) = tau * std::sqrt(yabs(t) * yabs(t) * y(t - 1) * y(t - 1));
    }
  } else {
    arma::vec x1 = yabs.subvec(1, Tsize - 1) - mua;
    arma::vec x2 = y.subvec(0, Tsize - 2) - muu;
    EH = arma::dot(x1, x2) / (Tsize - 2.0);
    for (int t = 1; t < Tsize; t++) {
      g_t(t, n_mom - 1) = yabs(t) * y(t - 1);
    }
  }

  // Aggregated moments
  double m1 = mu + 1.2704 - std::log(sigy * sigy);
  double m2 = gam0 + gamk_all(1) - (M_PI * M_PI / 2.0) -
    (1.0 / (1.0 - phi(0))) * (gamtmp + sigv * sigv);
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m4 = rho - (EH * std::sqrt(2.0 * M_PI)) /
    (sigv * sigy * sigy) * std::exp(-0.25 * gammatilde);

  // Per-obs scores
  g_t.col(0) = ly2 + 1.2704 - std::log(sigy * sigy);
  g_t.col(1) = gam0_t + gamk_t.col(0) - (M_PI * M_PI / 2.0) -
    (1.0 / (1.0 - phi(0))) * (gamtmp_t + sigv * sigv);

  arma::mat mk_t_mat = build_mk_t(phi, gamk_t, p);
  g_t.cols(2, p + 1) = mk_t_mat;

  // Leverage per-obs: rho - (g_t[,last] * sqrt(2pi)) / (sigv * sigy^2) * exp(-gt/4)
  double lev_scale = std::sqrt(2.0 * M_PI) / (sigv * sigy * sigy) *
    std::exp(-0.25 * gammatilde);
  for (int t = 1; t < Tsize; t++) {
    g_t(t, n_mom - 1) = rho - g_t(t, n_mom - 1) * lev_scale;
  }

  // Trim initial p rows
  g_t = g_t.rows(p, Tsize - 1);

  arma::vec g(n_mom);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m4;

  // HAC — BUG FIX: use T_eff and T_eff-1 consistently (was Tsize-3)
  Rcpp::Function pinv_f(pinv_fn.get());
  arma::mat Amat = hac_estimate(g_t, Tsize, Bartlett, pinv_f);

  return arma::as_scalar(g.t() * Amat * g);
}


// [[Rcpp::export]]
double LRT_moment_ar_Amat_cpp(arma::vec y, Rcpp::List mdl_out,
                                double del = 1e-10, bool Bartlett = true,
                                Rcpp::Nullable<Rcpp::Function> pinv_fn = R_NilValue) {
  int Tsize = y.n_elem;
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  int n_mom = p + 2;

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);
  double gamtmp = compute_gamtmp(phi, gamk_all.subvec(1, 2 * p), p);

  arma::mat gamk_t = build_gamk_t(ys, 2 * p);
  arma::vec gam0_t = ys % ys;
  arma::vec gamtmp_t = compute_gamtmp_t(phi, gamk_t, p);

  double m1 = mu + 1.2704 - std::log(sigy * sigy);
  double m2 = gam0 + gamk_all(1) - (M_PI * M_PI / 2.0) -
    (1.0 / (1.0 - phi(0))) * (gamtmp + sigv * sigv);
  arma::vec mk = build_mk(phi, gamk_all, p);

  arma::mat g_t(Tsize, n_mom, arma::fill::zeros);
  g_t.col(0) = ly2 + 1.2704 - std::log(sigy * sigy);
  g_t.col(1) = gam0_t + gamk_t.col(0) - (M_PI * M_PI / 2.0) -
    (1.0 / (1.0 - phi(0))) * (gamtmp_t + sigv * sigv);

  arma::mat mk_t_mat = build_mk_t(phi, gamk_t, p);
  g_t.cols(2, p + 1) = mk_t_mat;

  g_t = g_t.rows(p, Tsize - 1);

  arma::vec g(n_mom);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;

  Rcpp::Function pinv_f(pinv_fn.get());
  arma::mat Amat = hac_estimate(g_t, Tsize, Bartlett, pinv_f);

  return arma::as_scalar(g.t() * Amat * g);
}


// [[Rcpp::export]]
double LRT_moment_lev_t_Amat_cpp(arma::vec y, Rcpp::List mdl_out,
                                   std::string rho_type,
                                   double del = 1e-10, bool Bartlett = true,
                                   Rcpp::Nullable<Rcpp::Function> pinv_fn = R_NilValue) {
  int Tsize = y.n_elem;
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double rho = mdl_out["rho"];
  double gammatilde = mdl_out["gammatilde"];
  double nu = mdl_out["v"];
  int n_mom = p + 4;

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);

  arma::mat gamk_t = build_gamk_t(ys, 2 * p);
  arma::vec gam0_t = ys % ys;

  // EH
  arma::vec yabs = arma::abs(y);
  double muu = arma::mean(y.subvec(0, Tsize - 2));
  double mua = arma::mean(yabs.subvec(1, Tsize - 1));
  double EH;
  if (rho_type == "kendall") {
    arma::vec x1 = yabs.subvec(1, Tsize - 1) - mua;
    arma::vec x2 = y.subvec(0, Tsize - 2) - muu;
    double tau = kendall_corr(x1, x2);
    double var1 = arma::as_scalar(arma::var(yabs.subvec(1, Tsize - 1)));
    double var2 = arma::as_scalar(arma::var(y.subvec(0, Tsize - 2)));
    EH = tau * std::sqrt(var1 * var2);
  } else {
    arma::vec x1 = yabs.subvec(1, Tsize - 1) - mua;
    arma::vec x2 = y.subvec(0, Tsize - 2) - muu;
    EH = arma::dot(x1, x2) / (Tsize - 2.0);
  }

  Rcpp::Function ARMAacf("ARMAacf", Rcpp::Environment::namespace_env("stats"));
  arma::vec acf_vals = Rcpp::as<arma::vec>(ARMAacf(Rcpp::Named("ar") = phi,
                                                     Rcpp::Named("lag.max") = 1));
  double rho_w1 = acf_vals(1);

  double mu_log_sq_t = R::psigamma(0.5, 0) - R::psigamma(nu / 2.0, 0) + std::log(nu);
  double var_log_sq_t = R::psigamma(0.5, 1) + R::psigamma(nu / 2.0, 1);
  double Ct = std::exp(std::log(nu / 2.0) +
    2.0 * (R::lgammafn((nu - 1.0) / 2.0) - R::lgammafn(nu / 2.0)));

  // Aggregated moments
  double m1 = mu - mu_log_sq_t - std::log(sigy * sigy);
  double m2 = gam0 - arma::dot(phi, gamk_all.subvec(1, p)) - var_log_sq_t - sigv * sigv;
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m_prof = var_log_sq_t - gam0 + gamk_all(1) / rho_w1;
  double m_lev = rho - (std::sqrt(2.0 * M_PI) * EH) /
    (sigv * sigy * sigy * Ct) * std::exp(-0.25 * gammatilde);

  arma::vec g(n_mom);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m_prof;
  g(p + 3) = m_lev;

  // Per-obs scores
  arma::mat g_t(Tsize, n_mom, arma::fill::zeros);
  g_t.col(0) = ly2 - mu_log_sq_t - std::log(sigy * sigy);
  arma::vec g2 = gam0_t;
  for (int j = 0; j < p; j++) {
    g2 -= phi(j) * gamk_t.col(j);
  }
  g_t.col(1) = g2 - var_log_sq_t - sigv * sigv;
  arma::mat mk_t_mat = build_mk_t(phi, gamk_t, p);
  g_t.cols(2, p + 1) = mk_t_mat;
  g_t.col(p + 2) = var_log_sq_t - gam0_t + gamk_t.col(0) / rho_w1;

  // Leverage per-obs
  double lev_scale = std::sqrt(2.0 * M_PI) / (sigv * sigy * sigy * Ct) *
    std::exp(-0.25 * gammatilde);
  for (int t = 1; t < Tsize; t++) {
    g_t(t, n_mom - 1) = rho - yabs(t) * y(t - 1) * lev_scale;
  }

  g_t = g_t.rows(p, Tsize - 1);

  Rcpp::Function pinv_f(pinv_fn.get());
  arma::mat Amat = hac_estimate(g_t, Tsize, Bartlett, pinv_f);

  return arma::as_scalar(g.t() * Amat * g);
}


// [[Rcpp::export]]
double LRT_moment_lev_ged_Amat_cpp(arma::vec y, Rcpp::List mdl_out,
                                     std::string rho_type,
                                     double del, bool Bartlett,
                                     Rcpp::Nullable<Rcpp::Function> pinv_fn,
                                     arma::vec gh_nodes,
                                     arma::vec gh_weights) {
  int Tsize = y.n_elem;
  arma::vec ly2 = arma::log(arma::square(y) + del);
  double mu = arma::mean(ly2);
  arma::vec ys = ly2 - mu;
  arma::vec phi = Rcpp::as<arma::vec>(mdl_out["phi"]);
  int p = phi.n_elem;
  double sigy = mdl_out["sigy"];
  double sigv = mdl_out["sigv"];
  double rho = mdl_out["rho"];
  double gammatilde = mdl_out["gammatilde"];
  double nu = mdl_out["v"];
  int n_mom = p + 4;

  arma::vec gamk_all = batch_acov(ys, 2 * p);
  double gam0 = gamk_all(0);

  arma::mat gamk_t = build_gamk_t(ys, 2 * p);
  arma::vec gam0_t = ys % ys;

  // EH
  arma::vec yabs = arma::abs(y);
  double muu = arma::mean(y.subvec(0, Tsize - 2));
  double mua = arma::mean(yabs.subvec(1, Tsize - 1));
  double EH;
  if (rho_type == "kendall") {
    arma::vec x1 = yabs.subvec(1, Tsize - 1) - mua;
    arma::vec x2 = y.subvec(0, Tsize - 2) - muu;
    double tau = kendall_corr(x1, x2);
    double var1 = arma::as_scalar(arma::var(yabs.subvec(1, Tsize - 1)));
    double var2 = arma::as_scalar(arma::var(y.subvec(0, Tsize - 2)));
    EH = tau * std::sqrt(var1 * var2);
  } else {
    arma::vec x1 = yabs.subvec(1, Tsize - 1) - mua;
    arma::vec x2 = y.subvec(0, Tsize - 2) - muu;
    EH = arma::dot(x1, x2) / (Tsize - 2.0);
  }

  Rcpp::Function ARMAacf("ARMAacf", Rcpp::Environment::namespace_env("stats"));
  arma::vec acf_vals = Rcpp::as<arma::vec>(ARMAacf(Rcpp::Named("ar") = phi,
                                                     Rcpp::Named("lag.max") = 1));
  double rho_w1 = acf_vals(1);

  double mu_log_sq_ged = (2.0 / nu) * R::psigamma(1.0 / nu, 0) +
    R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu);
  double var_log_sq_ged = std::pow(2.0 / nu, 2) * R::psigamma(1.0 / nu, 1);
  double E_abs_u = std::exp(0.5 * (R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu)) +
    R::lgammafn(2.0 / nu) - R::lgammafn(1.0 / nu));

  double c_shift = rho * sigv / 2.0;
  double E_g_shifted = 0.0;
  for (int i = 0; i < (int)gh_nodes.n_elem; i++) {
    double p_val = R::pnorm5(gh_nodes(i) + c_shift, 0.0, 1.0, 1, 0);
    E_g_shifted += gh_weights(i) * qged_std_cpp(p_val, nu);
  }

  // Aggregated
  double m1 = mu - mu_log_sq_ged - std::log(sigy * sigy);
  double m2 = gam0 - arma::dot(phi, gamk_all.subvec(1, p)) - var_log_sq_ged - sigv * sigv;
  arma::vec mk = build_mk(phi, gamk_all, p);
  double m_prof = var_log_sq_ged - gam0 + gamk_all(1) / rho_w1;
  double m_lev = sigy * sigy * E_abs_u * E_g_shifted *
    std::exp(gammatilde / 4.0) - EH;

  arma::vec g(n_mom);
  g(0) = m1;
  g(1) = m2;
  g.subvec(2, p + 1) = mk;
  g(p + 2) = m_prof;
  g(p + 3) = m_lev;

  // Per-obs scores
  arma::mat g_t(Tsize, n_mom, arma::fill::zeros);
  g_t.col(0) = ly2 - mu_log_sq_ged - std::log(sigy * sigy);
  arma::vec g2 = gam0_t;
  for (int j = 0; j < p; j++) {
    g2 -= phi(j) * gamk_t.col(j);
  }
  g_t.col(1) = g2 - var_log_sq_ged - sigv * sigv;
  arma::mat mk_t_mat = build_mk_t(phi, gamk_t, p);
  g_t.cols(2, p + 1) = mk_t_mat;
  g_t.col(p + 2) = var_log_sq_ged - gam0_t + gamk_t.col(0) / rho_w1;

  // Leverage per-obs
  double lev_const = sigy * sigy * E_abs_u * E_g_shifted * std::exp(gammatilde / 4.0);
  for (int t = 1; t < Tsize; t++) {
    g_t(t, n_mom - 1) = lev_const - yabs(t) * y(t - 1);
  }

  g_t = g_t.rows(p, Tsize - 1);

  Rcpp::Function pinv_f(pinv_fn.get());
  arma::mat Amat = hac_estimate(g_t, Tsize, Bartlett, pinv_f);

  return arma::as_scalar(g.t() * Amat * g);
}
