// ============================================================================ //
// wARMASVp: Kalman Filter and GMKF for SV(p) Models
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declaration: GED CDF, defined in utils_cpp.cpp
double pged_std_cpp(double u, double nu);

// --------------------------------------------------------------------------- //
// Compute the leverage-shift proxy zeta_hat from the observed u, by error type.
// dist_code:    0 = Gaussian, 1 = Student-t, 2 = GED
// proxy_type:   0 = "u" (paper-faithful per SVHT Remark 3.5)
//               1 = "bayes_optimal" (Bayes posterior mean for Student-t)
// For Gaussian, u = zeta exactly (no choice).
// For GED, u = qged_std(Phi(zeta)) so zeta = Phi^{-1}(F_GED(u;nu)) (exact;
//   the proxy_type flag has no effect).
// For Student-t, u = zeta * lambda^{-1/2}.  The "u" proxy uses zeta_hat = u
//   (paper choice).  The "bayes_optimal" proxy uses
//     E[zeta | u] = u * sqrt((nu+1)/(nu+u^2)) * Gamma((nu+2)/2)/Gamma((nu+1)/2)
//                * sqrt(2/(nu+1))
//             = u * sqrt(2/(nu+u^2)) * Gamma((nu+2)/2)/Gamma((nu+1)/2).
// --------------------------------------------------------------------------- //
static inline double leverage_zeta_proxy(double u, int dist_code, double nu,
                                          int proxy_type) {
  if (dist_code == 0) {
    return u;                       // Gaussian: zeta = u exactly
  } else if (dist_code == 2) {
    double p = pged_std_cpp(u, nu); // GED: invert the copula
    return R::qnorm(p, 0.0, 1.0, 1, 0);
  } else {
    // Student-t
    if (proxy_type == 0) {
      return u;                     // paper-faithful (Remark 3.5)
    } else {
      // Bayes-optimal: posterior mean of zeta given u
      double factor = std::sqrt(2.0 / (nu + u * u))
                    * std::exp(R::lgammafn((nu + 2.0) / 2.0)
                             - R::lgammafn((nu + 1.0) / 2.0));
      return u * factor;
    }
  }
}

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
                       const arma::mat& P0,
                       int dist_code = 0,
                       double nu = 0.0,
                       int proxy_type = 0) {
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
        // Apply distribution-aware leverage proxy: zeta_hat = u for Gaussian
        // (exact), Phi^{-1}(F_GED(u;nu)) for GED (exact), and either u
        // (paper-faithful) or E[zeta|u] (Bayes-optimal) for Student-t.
        double zeta_hat = leverage_zeta_proxy(zt(t - 1), dist_code, nu,
                                               proxy_type);
        xi_p += sigma_v * delta_p * zeta_hat * h_vec;
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
    arma::mat A = P_filt_arr.slice(t) * F_mat.t() * arma::pinv(P_pred_arr.slice(t + 1));
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
    Named("P_filt_T") = P_filt_arr.slice(T - 1),
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
                     const arma::mat& P0,
                     int dist_code = 0,
                     double nu = 0.0,
                     int proxy_type = 0) {
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
        double zeta_hat = leverage_zeta_proxy(zt(t - 1), dist_code, nu,
                                               proxy_type);
        xi_p += sigma_v * delta_p * zeta_hat * h_vec;
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
    arma::mat A = P_filt_arr.slice(t) * F_mat.t() * arma::pinv(P_pred_arr.slice(t + 1));
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
    Named("P_filt_T") = P_filt_arr.slice(T - 1),
    Named("xi_filtered") = xi_filt,
    Named("xi_smoothed") = xi_smooth,
    Named("loglik") = loglik
  );
}


// --------------------------------------------------------------------------- //
// EM loop for Gaussian-mixture fit to measurement-noise density (KSC 1998).
// R-side supplies initial (m, s2, q) from quantile binning; this function
// iterates E/M until convergence and returns sorted mixture parameters.
// Pointer-arithmetic hot loops in E-step and M-step for ~15-20x over the R
// implementation at n = 1e5 (dominated by row-max in log-sum-exp normalize).
// --------------------------------------------------------------------------- //
// [[Rcpp::export]]
List fit_ksc_em_cpp(const arma::vec& eps,
                    const arma::vec& init_m,
                    const arma::vec& init_s2,
                    const arma::vec& init_q,
                    int max_iter = 500,
                    double tol = 1e-8) {
  const arma::uword n = eps.n_elem;
  const arma::uword K = init_m.n_elem;
  arma::vec m = init_m;
  arma::vec s2 = init_s2;
  arma::vec q = init_q;
  arma::mat log_tau(n, K);
  arma::vec max_log(n);
  arma::mat tau(n, K);
  arma::vec row_sums(n);
  arma::vec q_new(K), m_new(K), s2_new(K);

  const double LOG_2PI = std::log(2.0 * arma::datum::pi);
  int n_iter_used = max_iter;

  for (int iter = 1; iter <= max_iter; ++iter) {
    // E-step: log_tau[i,k] = log q_k - 0.5 log(2 pi s2_k) - 0.5 (eps_i - m_k)^2 / s2_k
    for (arma::uword k = 0; k < K; ++k) {
      const double inv_s2 = 1.0 / s2(k);
      const double c = std::log(q(k)) - 0.5 * (LOG_2PI + std::log(s2(k)));
      const double mk = m(k);
      double* p_out = log_tau.colptr(k);
      const double* p_eps = eps.memptr();
      for (arma::uword i = 0; i < n; ++i) {
        const double d = p_eps[i] - mk;
        p_out[i] = c - 0.5 * d * d * inv_s2;
      }
    }
    // Normalise row-wise via log-sum-exp
    for (arma::uword i = 0; i < n; ++i) {
      double mx = log_tau(i, 0);
      for (arma::uword k = 1; k < K; ++k)
        if (log_tau(i, k) > mx) mx = log_tau(i, k);
      max_log(i) = mx;
    }
    for (arma::uword k = 0; k < K; ++k) {
      double* p_tau = tau.colptr(k);
      const double* p_lt = log_tau.colptr(k);
      for (arma::uword i = 0; i < n; ++i)
        p_tau[i] = std::exp(p_lt[i] - max_log(i));
    }
    row_sums.zeros();
    for (arma::uword k = 0; k < K; ++k) {
      double* p_rs = row_sums.memptr();
      const double* p_tau = tau.colptr(k);
      for (arma::uword i = 0; i < n; ++i) p_rs[i] += p_tau[i];
    }
    for (arma::uword k = 0; k < K; ++k) {
      double* p_tau = tau.colptr(k);
      const double* p_rs = row_sums.memptr();
      for (arma::uword i = 0; i < n; ++i) p_tau[i] /= p_rs[i];
    }

    // M-step
    for (arma::uword k = 0; k < K; ++k) {
      double Nk = 0.0, s_eps = 0.0;
      const double* p_tau = tau.colptr(k);
      const double* p_eps = eps.memptr();
      for (arma::uword i = 0; i < n; ++i) {
        Nk    += p_tau[i];
        s_eps += p_tau[i] * p_eps[i];
      }
      if (Nk < 1e-10) Nk = 1e-10;
      q_new(k) = Nk / (double)n;
      m_new(k) = s_eps / Nk;
      double s_var = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        const double d = p_eps[i] - m_new(k);
        s_var += p_tau[i] * d * d;
      }
      s2_new(k) = s_var / Nk;
      if (s2_new(k) < 1e-10) s2_new(k) = 1e-10;
    }

    const double dm = arma::abs(m_new - m).max();
    const double dq = arma::abs(q_new - q).max();
    q = q_new; m = m_new; s2 = s2_new;
    if (dm < tol && dq < tol) { n_iter_used = iter; break; }
  }
  arma::uvec ord = arma::sort_index(m);
  arma::vec q_out = q(ord);
  arma::vec m_out = m(ord);
  arma::vec s_out = s2(ord);
  return List::create(Named("weights") = q_out,
                      Named("means")   = m_out,
                      Named("vars")    = s_out,
                      Named("n_iter")  = n_iter_used);
}
