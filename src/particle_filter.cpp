// ============================================================================ //
// wARMASVp: Bootstrap Particle Filter for SV(p) Models
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declarations
arma::mat solve_lyapunov_discrete_cpp(const arma::mat& F_mat, const arma::mat& Q);
double qged_std_cpp(double p, double nu);

// --------------------------------------------------------------------------- //
// Exact measurement noise log-densities (centered)
// --------------------------------------------------------------------------- //

// Log-density of centered log-chi-squared(1) (Gaussian case)
// eps = log(z^2) - E[log(z^2)] where z ~ N(0,1)
static double log_density_eps_gaussian(double eps) {
  double mu_bar = R::digamma(0.5) + std::log(2.0);
  double x = eps + mu_bar;  // raw log(chi2_1) scale
  // log(chi2_1) density: f(x) = (1/sqrt(2*pi)) * exp(x/2 - exp(x)/2)
  return -0.5 * std::log(2.0 * M_PI) + x / 2.0 - std::exp(x) / 2.0;
}

// Log-density of centered log-F(1,nu) (Student-t case)
static double log_density_eps_t(double eps, double nu, double mu_bar) {
  double x = eps + mu_bar;
  return -0.5 * std::log(nu) - R::lbeta(0.5, nu / 2.0) +
         x / 2.0 - ((1.0 + nu) / 2.0) * std::log(1.0 + std::exp(x) / nu);
}

// Log-density of centered log-GED^2 (GED case)
static double log_density_eps_ged(double eps, double nu, double a, double mu_bar_g) {
  double x = eps + mu_bar_g;
  double u_abs = std::exp(x / 2.0);
  return std::log(nu) - std::log(2.0) - std::log(a) - R::lgammafn(1.0 / nu) -
         std::pow(u_abs / a, nu) + x / 2.0;
}

// pged_std_cpp moved to utils_cpp.cpp (now shared with kalman_filter.cpp).
// Forward declaration only.
double pged_std_cpp(double u, double nu);

// --------------------------------------------------------------------------- //
// Systematic resampling
// --------------------------------------------------------------------------- //
static arma::uvec systematic_resample(const arma::vec& weights, int M) {
  arma::uvec indices(M);
  double u0 = R::runif(0.0, 1.0 / M);
  double cumsum = weights(0);
  int j = 0;
  for (int i = 0; i < M; i++) {
    double target = u0 + (double)i / (double)M;
    while (cumsum < target && j < M - 1) {
      j++;
      cumsum += weights(j);
    }
    indices(i) = j;
  }
  return indices;
}

// --------------------------------------------------------------------------- //
// Bootstrap Particle Filter
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
List particle_filter_svp_cpp(
    const arma::vec& y_raw,
    const arma::vec& phi,
    double sigma_y,
    double sigma_v,
    double nu,
    int dist_code,       // 0 = gaussian, 1 = student_t, 2 = ged
    double delta,
    int M,
    int seed,
    double del) {

  int T_obs = y_raw.n_elem;
  int p = phi.n_elem;

  // Set R's RNG seed
  Environment base("package:base");
  Function set_seed = base["set.seed"];
  set_seed(seed);

  // Build companion matrix F
  arma::mat F_mat(p, p, arma::fill::zeros);
  F_mat.row(0) = phi.t();
  if (p > 1) {
    for (int j = 0; j < p - 1; j++) {
      F_mat(j + 1, j) = 1.0;
    }
  }
  arma::vec h_vec(p, arma::fill::zeros);
  h_vec(0) = 1.0;
  arma::vec r_vec = h_vec;

  // Pre-compute distribution constants
  double mu_bar = 0.0, a_ged = 0.0, mu_bar_g = 0.0;
  if (dist_code == 0) {
    mu_bar = R::digamma(0.5) + std::log(2.0);
  } else if (dist_code == 1) {
    mu_bar = R::digamma(0.5) - R::digamma(nu / 2.0) + std::log(nu);
  } else {
    a_ged = std::exp(0.5 * (R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu)));
    mu_bar_g = (2.0 / nu) * R::digamma(1.0 / nu) +
               R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu);
  }
  double mu_intercept = std::log(sigma_y * sigma_y) + mu_bar;
  if (dist_code == 2) {
    mu_intercept = std::log(sigma_y * sigma_y) + mu_bar_g;
  }

  // Initial state covariance (Lyapunov)
  arma::mat Q_init = sigma_v * sigma_v * (r_vec * r_vec.t());
  arma::mat P0 = solve_lyapunov_discrete_cpp(F_mat, Q_init);

  // Draw initial particles from N(0, P0) via Cholesky
  arma::mat L0 = arma::chol(P0, "lower");
  arma::mat particles(p, M);  // p x M
  for (int i = 0; i < M; i++) {
    arma::vec z = Rcpp::as<arma::vec>(Rcpp::rnorm(p, 0.0, 1.0));
    particles.col(i) = L0 * z;
  }

  // Storage for z_t particles (leverage)
  arma::vec z_particles(M, arma::fill::zeros);

  // Leverage parameters
  double delta_eff = (delta != 0.0) ? std::sqrt(1.0 - delta * delta) : 0.0;

  // Output storage
  arma::vec w_filt(T_obs), P_filt(T_obs), ess_vec(T_obs);
  double loglik = 0.0;

  // Full p x p filtered state and covariance at final time step (for forecasting)
  arma::vec xi_filt_T(p, arma::fill::zeros);
  arma::mat P_filt_T_mat(p, p, arma::fill::zeros);

  // Log-squared observations
  arma::vec y_star(T_obs);
  for (int t = 0; t < T_obs; t++) {
    y_star(t) = std::log(y_raw(t) * y_raw(t) + del);
  }

  // --- Main filter loop ---
  for (int t = 0; t < T_obs; t++) {

    // PROPAGATE
    for (int i = 0; i < M; i++) {
      double v_t;
      if (delta == 0.0) {
        v_t = R::rnorm(0.0, sigma_v);
      } else {
        double eta = R::rnorm(0.0, 1.0);
        v_t = sigma_v * (delta * z_particles(i) + delta_eff * eta);
      }
      arma::vec xi_new = F_mat * particles.col(i) + v_t * r_vec;
      particles.col(i) = xi_new;
    }

    // WEIGHT
    arma::vec log_weights(M);
    for (int i = 0; i < M; i++) {
      double w_i = arma::dot(h_vec, particles.col(i));
      double eps_i = y_star(t) - mu_intercept - w_i;

      if (dist_code == 0) {
        log_weights(i) = log_density_eps_gaussian(eps_i);
      } else if (dist_code == 1) {
        log_weights(i) = log_density_eps_t(eps_i, nu, mu_bar);
      } else {
        log_weights(i) = log_density_eps_ged(eps_i, nu, a_ged, mu_bar_g);
      }
    }

    // Normalize (log-sum-exp)
    double max_lw = log_weights.max();
    arma::vec weights = arma::exp(log_weights - max_lw);
    double sum_w = arma::sum(weights);
    weights /= sum_w;
    loglik += max_lw + std::log(sum_w) - std::log((double)M);

    // ESS (before resampling)
    ess_vec(t) = 1.0 / arma::dot(weights, weights);

    // Weighted mean and variance (before resampling)
    arma::vec w_vals(M);
    for (int i = 0; i < M; i++) {
      w_vals(i) = arma::dot(h_vec, particles.col(i));
    }
    double w_mean = arma::dot(weights, w_vals);
    double w_var = 0.0;
    for (int i = 0; i < M; i++) {
      double d = w_vals(i) - w_mean;
      w_var += weights(i) * d * d;
    }
    w_filt(t) = w_mean;
    P_filt(t) = w_var;

    // At final time step, compute full p x p filtered state and covariance
    if (t == T_obs - 1) {
      xi_filt_T.zeros();
      P_filt_T_mat.zeros();
      for (int i = 0; i < M; i++) {
        xi_filt_T += weights(i) * particles.col(i);
      }
      for (int i = 0; i < M; i++) {
        arma::vec diff = particles.col(i) - xi_filt_T;
        P_filt_T_mat += weights(i) * (diff * diff.t());
      }
    }

    // RESAMPLE (systematic)
    arma::uvec idx = systematic_resample(weights, M);
    arma::mat new_particles(p, M);
    for (int i = 0; i < M; i++) {
      new_particles.col(i) = particles.col(idx(i));
    }
    particles = new_particles;

    // RECOVER z_t (for leverage at t+1)
    if (delta != 0.0) {
      for (int i = 0; i < M; i++) {
        double w_i = arma::dot(h_vec, particles.col(i));
        double u_i = y_raw(t) / (sigma_y * std::exp(w_i / 2.0));

        if (dist_code == 0) {
          // Gaussian: u = z (exact)
          z_particles(i) = u_i;
        } else if (dist_code == 1) {
          // Student-t: sample lambda from POSTERIOR given u_i.  Under the
          // scale-mixture u_i = zeta_i * lambda_i^{-1/2} with zeta ~ N(0,1)
          // and lambda ~ Gamma(nu/2, rate=nu/2), Bayes' rule gives
          //   lambda | u_i ~ Gamma((nu+1)/2, rate=(nu + u_i^2)/2).
          // Then zeta = u_i * sqrt(lambda) is the conditional draw.
          // (Previous code sampled from the prior Gamma(nu/2, nu/2), which
          // overshoots leverage feedback in the tails.  Bug fix 2026-05-09.)
          double lambda = R::rgamma((nu + 1.0) / 2.0,
                                     2.0 / (nu + u_i * u_i));
          z_particles(i) = u_i * std::sqrt(lambda);
        } else {
          // GED: z = Phi^{-1}(F_GED(u))
          double p_ged = pged_std_cpp(u_i, nu);
          // Clamp to avoid Inf
          if (p_ged < 1e-10) p_ged = 1e-10;
          if (p_ged > 1.0 - 1e-10) p_ged = 1.0 - 1e-10;
          z_particles(i) = R::qnorm(p_ged, 0.0, 1.0, 1, 0);
        }
      }
    }
  }

  return List::create(
    Named("w_filtered") = w_filt,
    Named("P_filtered") = P_filt,
    Named("ESS") = ess_vec,
    Named("loglik") = loglik,
    Named("xi_filt_T") = xi_filt_T,
    Named("P_filt_T") = P_filt_T_mat
  );
}
