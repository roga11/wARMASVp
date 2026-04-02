// ============================================================================ //
// wARMASVp: Kalman Filter and GMKF for SV(p) Models
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// --------------------------------------------------------------------------- //
// Discrete Lyapunov solver: X = F X F' + Q
// Uses vectorization: vec(X) = (I - F ⊗ F)^{-1} vec(Q)
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
arma::mat solve_lyapunov_discrete_cpp(const arma::mat& F_mat, const arma::mat& Q) {
  int p = F_mat.n_rows;
  arma::mat I_p2 = arma::eye(p * p, p * p);
  arma::mat FkF = arma::kron(F_mat, F_mat);
  arma::vec Q_vec = arma::vectorise(Q);
  arma::vec X_vec = arma::solve(I_p2 - FkF, Q_vec);
  return arma::reshape(X_vec, p, p);
}

// --------------------------------------------------------------------------- //
// Corrected Kalman Filter (CKF) — forward pass + RTS smoother
//
// This is the standard Kalman filter using distribution-specific σ_ε²(ν)
// as measurement noise variance, with leverage support.
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
List kalman_filter_cpp(const arma::vec& y_star,
                       const arma::vec& y_raw,
                       const arma::vec& phi,
                       double sigma_y,
                       double sigma_v,
                       double delta_p,
                       double sig_eps2,
                       double var_zt,
                       const arma::mat& P0) {
  int p = phi.n_elem;
  int T = y_star.n_elem;

  // Build companion matrix F
  arma::mat F_mat(p, p, arma::fill::zeros);
  F_mat.row(0) = phi.t();
  if (p > 1) {
    for (int j = 0; j < p - 1; j++) {
      F_mat(j + 1, j) = 1.0;
    }
  }

  // Selection vectors
  arma::vec h_vec(p, arma::fill::zeros);
  h_vec(0) = 1.0;

  // State noise covariance (no leverage component, handled in prediction)
  double sigma_v_eff2 = (delta_p == 0.0) ? sigma_v * sigma_v
                                          : sigma_v * sigma_v * (1.0 - delta_p * delta_p);
  double sigma_v_lev2 = sigma_v * sigma_v * delta_p * delta_p * var_zt;
  arma::mat Q_eff = sigma_v_eff2 * (h_vec * h_vec.t());
  arma::mat Q_lev = sigma_v_lev2 * (h_vec * h_vec.t());

  // Storage
  arma::mat xi_pred(p, T, arma::fill::zeros);
  arma::mat xi_filt(p, T, arma::fill::zeros);
  arma::cube P_pred_arr(p, p, T, arma::fill::zeros);
  arma::cube P_filt_arr(p, p, T, arma::fill::zeros);
  arma::vec w_filt(T);
  arma::vec w_pred(T);
  arma::vec P_filt_11(T);
  arma::vec P_pred_11(T);
  arma::vec zt(T);
  double loglik = 0.0;

  // Initialize
  arma::vec xi_hat = arma::zeros(p);
  arma::mat P_hat = P0;

  // Forward pass
  for (int t = 0; t < T; t++) {
    // Prediction
    arma::vec xi_p;
    arma::mat P_p;
    if (t == 0) {
      xi_p = F_mat * xi_hat;
      P_p = F_mat * P_hat * F_mat.t() + sigma_v * sigma_v * (h_vec * h_vec.t());
    } else {
      xi_p = F_mat * xi_filt.col(t - 1);
      if (delta_p != 0.0) {
        xi_p += sigma_v * delta_p * zt(t - 1) * h_vec;
      }
      P_p = F_mat * P_filt_arr.slice(t - 1) * F_mat.t() + Q_eff + Q_lev;
    }
    xi_pred.col(t) = xi_p;
    P_pred_arr.slice(t) = P_p;
    w_pred(t) = arma::dot(h_vec, xi_p);
    P_pred_11(t) = P_p(0, 0);

    // Update
    double innov = y_star(t) - arma::dot(h_vec, xi_p);
    double f_t = arma::as_scalar(h_vec.t() * P_p * h_vec) + sig_eps2;
    arma::vec K_gain = P_p * h_vec / f_t;

    arma::vec xi_u = xi_p + K_gain * innov;
    arma::mat P_u = P_p - K_gain * h_vec.t() * P_p;

    xi_filt.col(t) = xi_u;
    P_filt_arr.slice(t) = P_u;
    w_filt(t) = arma::dot(h_vec, xi_u);
    P_filt_11(t) = P_u(0, 0);

    // Standardized residual
    zt(t) = (y_raw(t) / sigma_y) * std::exp(-0.5 * w_filt(t));

    // Log-likelihood contribution
    loglik += -0.5 * (std::log(2.0 * M_PI * f_t) + innov * innov / f_t);
  }

  // RTS Backward pass (smoothing)
  arma::mat xi_smooth(p, T, arma::fill::zeros);
  arma::cube P_smooth_arr(p, p, T, arma::fill::zeros);
  xi_smooth.col(T - 1) = xi_filt.col(T - 1);
  P_smooth_arr.slice(T - 1) = P_filt_arr.slice(T - 1);

  for (int t = T - 2; t >= 0; t--) {
    arma::mat A = P_filt_arr.slice(t) * F_mat.t() * arma::inv(P_pred_arr.slice(t + 1));
    xi_smooth.col(t) = xi_filt.col(t) + A * (xi_smooth.col(t + 1) - xi_pred.col(t + 1));
    P_smooth_arr.slice(t) = P_filt_arr.slice(t) +
      A * (P_smooth_arr.slice(t + 1) - P_pred_arr.slice(t + 1)) * A.t();
  }

  // Extract smoothed w and zt_smoothed
  arma::vec w_smooth(T);
  arma::vec zt_smooth(T);
  for (int t = 0; t < T; t++) {
    w_smooth(t) = arma::dot(h_vec, xi_smooth.col(t));
    zt_smooth(t) = (y_raw(t) / sigma_y) * std::exp(-0.5 * w_smooth(t));
  }

  return List::create(
    Named("w_filtered") = w_filt,
    Named("w_smoothed") = w_smooth,
    Named("w_predicted") = w_pred,
    Named("zt") = zt,
    Named("zt_smoothed") = zt_smooth,
    Named("P_filtered") = P_filt_11,
    Named("P_predicted") = P_pred_11,
    Named("xi_filtered") = xi_filt,
    Named("xi_smoothed") = xi_smooth,
    Named("loglik") = loglik
  );
}

// --------------------------------------------------------------------------- //
// Gaussian Mixture Kalman Filter (GMKF)
//
// Runs K parallel Kalman updates per time step (one per mixture component),
// then collapses via moment matching. Includes RTS smoothing on collapsed states.
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
List gmkf_filter_cpp(const arma::vec& y_star,
                     const arma::vec& y_raw,
                     const arma::vec& phi,
                     double sigma_y,
                     double sigma_v,
                     double delta_p,
                     double var_zt,
                     const arma::vec& mix_weights,
                     const arma::vec& mix_means,
                     const arma::vec& mix_vars,
                     double mu_intercept,
                     const arma::mat& P0) {
  int p = phi.n_elem;
  int T = y_star.n_elem;
  int K = mix_weights.n_elem;

  // Build companion matrix F
  arma::mat F_mat(p, p, arma::fill::zeros);
  F_mat.row(0) = phi.t();
  if (p > 1) {
    for (int j = 0; j < p - 1; j++) {
      F_mat(j + 1, j) = 1.0;
    }
  }

  // Selection vectors
  arma::vec h_vec(p, arma::fill::zeros);
  h_vec(0) = 1.0;

  // State noise covariance
  double sigma_v_eff2 = (delta_p == 0.0) ? sigma_v * sigma_v
                                          : sigma_v * sigma_v * (1.0 - delta_p * delta_p);
  double sigma_v_lev2 = sigma_v * sigma_v * delta_p * delta_p * var_zt;
  arma::mat Q_eff = sigma_v_eff2 * (h_vec * h_vec.t());
  arma::mat Q_lev = sigma_v_lev2 * (h_vec * h_vec.t());

  // Storage
  arma::mat xi_pred(p, T, arma::fill::zeros);
  arma::mat xi_filt(p, T, arma::fill::zeros);
  arma::cube P_pred_arr(p, p, T, arma::fill::zeros);
  arma::cube P_filt_arr(p, p, T, arma::fill::zeros);
  arma::vec w_filt(T), w_pred(T), P_filt_11(T), P_pred_11(T), zt(T);
  double loglik = 0.0;

  // Initialize
  arma::vec xi_hat = arma::zeros(p);
  arma::mat P_hat = P0;

  for (int t = 0; t < T; t++) {
    // --- Prediction ---
    arma::vec xi_p;
    arma::mat P_p;
    if (t == 0) {
      xi_p = F_mat * xi_hat;
      P_p = F_mat * P_hat * F_mat.t() + sigma_v * sigma_v * (h_vec * h_vec.t());
    } else {
      xi_p = F_mat * xi_filt.col(t - 1);
      if (delta_p != 0.0) {
        xi_p += sigma_v * delta_p * zt(t - 1) * h_vec;
      }
      P_p = F_mat * P_filt_arr.slice(t - 1) * F_mat.t() + Q_eff + Q_lev;
    }
    xi_pred.col(t) = xi_p;
    P_pred_arr.slice(t) = P_p;
    w_pred(t) = arma::dot(h_vec, xi_p);
    P_pred_11(t) = P_p(0, 0);

    // --- Per-component update ---
    arma::mat xi_k(p, K);
    arma::cube P_k(p, p, K);
    arma::vec log_w(K);

    for (int k = 0; k < K; k++) {
      // Innovation for component k: y*_t - (mu + m_k) - h' xi_p
      // Note: y_star already = log(y²+del) - mu, so the intercept is 0 for CKF
      // For GMKF: y_star = log(y²) raw, and mu_intercept + m_k is subtracted
      double e_tk = y_star(t) - mu_intercept - mix_means(k) - arma::dot(h_vec, xi_p);
      double f_tk = arma::as_scalar(h_vec.t() * P_p * h_vec) + mix_vars(k);
      arma::vec K_gain = P_p * h_vec / f_tk;

      xi_k.col(k) = xi_p + K_gain * e_tk;
      P_k.slice(k) = P_p - K_gain * h_vec.t() * P_p;

      log_w(k) = std::log(mix_weights(k)) - 0.5 * (std::log(f_tk) + e_tk * e_tk / f_tk);
    }

    // Normalize weights (log-sum-exp)
    double max_lw = log_w.max();
    arma::vec w_norm = arma::exp(log_w - max_lw);
    double sum_w = arma::sum(w_norm);
    w_norm /= sum_w;
    loglik += max_lw + std::log(sum_w);

    // --- Collapse (moment matching) ---
    arma::vec xi_u = xi_k * w_norm;
    arma::mat P_u(p, p, arma::fill::zeros);
    for (int k = 0; k < K; k++) {
      arma::vec diff_k = xi_k.col(k) - xi_u;
      P_u += w_norm(k) * (P_k.slice(k) + diff_k * diff_k.t());
    }

    xi_filt.col(t) = xi_u;
    P_filt_arr.slice(t) = P_u;
    w_filt(t) = arma::dot(h_vec, xi_u);
    P_filt_11(t) = P_u(0, 0);

    // Standardized residual
    zt(t) = (y_raw(t) / sigma_y) * std::exp(-0.5 * w_filt(t));
  }

  // --- RTS Backward pass (collapsed smoother, Kim 1994 approximation) ---
  arma::mat xi_smooth(p, T, arma::fill::zeros);
  arma::cube P_smooth_arr(p, p, T, arma::fill::zeros);
  xi_smooth.col(T - 1) = xi_filt.col(T - 1);
  P_smooth_arr.slice(T - 1) = P_filt_arr.slice(T - 1);

  for (int t = T - 2; t >= 0; t--) {
    arma::mat P_pred_inv = arma::inv(P_pred_arr.slice(t + 1));
    arma::mat A = P_filt_arr.slice(t) * F_mat.t() * P_pred_inv;
    xi_smooth.col(t) = xi_filt.col(t) + A * (xi_smooth.col(t + 1) - xi_pred.col(t + 1));
    P_smooth_arr.slice(t) = P_filt_arr.slice(t) +
      A * (P_smooth_arr.slice(t + 1) - P_pred_arr.slice(t + 1)) * A.t();
  }

  // Extract smoothed w and zt_smoothed
  arma::vec w_smooth(T);
  arma::vec zt_smooth(T);
  for (int t = 0; t < T; t++) {
    w_smooth(t) = arma::dot(h_vec, xi_smooth.col(t));
    zt_smooth(t) = (y_raw(t) / sigma_y) * std::exp(-0.5 * w_smooth(t));
  }

  return List::create(
    Named("w_filtered") = w_filt,
    Named("w_smoothed") = w_smooth,
    Named("w_predicted") = w_pred,
    Named("zt") = zt,
    Named("zt_smoothed") = zt_smooth,
    Named("P_filtered") = P_filt_11,
    Named("P_predicted") = P_pred_11,
    Named("xi_filtered") = xi_filt,
    Named("xi_smoothed") = xi_smooth,
    Named("loglik") = loglik
  );
}
