# =========================================================================== #
# Kalman Filtering, GMKF, and Filtering API for SV(p) Models
# =========================================================================== #


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

#' Solve Discrete Lyapunov Equation
#'
#' Solves X = F X t(F) + Q for X using the vectorization approach.
#'
#' @param F_mat Square matrix.
#' @param Q Square matrix (same dimensions as \code{F_mat}).
#' @return Solution matrix \code{X}.
#' @keywords internal
solve_lyapunov_discrete <- function(F_mat, Q) {
  solve_lyapunov_discrete_cpp(F_mat, Q)
}

#' CDF of Standardized GED
#'
#' Computes the CDF of the standardized GED(\eqn{\nu}) distribution with
#' unit variance.
#'
#' @param u Numeric vector. Evaluation points.
#' @param nu Numeric. GED shape parameter (\eqn{\nu > 0}).
#' @return Numeric vector of probabilities.
#' @keywords internal
pged_std <- function(u, nu) {
  a <- exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu)))
  ifelse(u >= 0,
         0.5 + 0.5 * pgamma((u / a)^nu, shape = 1 / nu),
         0.5 - 0.5 * pgamma((-u / a)^nu, shape = 1 / nu))
}


# --------------------------------------------------------------------------- #
# Measurement noise density functions (for verification and future BPF)
# --------------------------------------------------------------------------- #

#' Density of Centered log-F(1,nu) Measurement Noise (Student-t)
#'
#' @param y Numeric vector. Evaluation points (centered: E[eps] = 0).
#' @param nu Numeric. Student-t degrees of freedom.
#' @return Density values.
#' @keywords internal
density_eps_t <- function(y, nu) {
  mu_bar <- digamma(0.5) - digamma(nu / 2) + log(nu)
  x <- y + mu_bar
  log_dens <- -0.5 * log(nu) - lbeta(0.5, nu / 2) +
    x / 2 - ((1 + nu) / 2) * log(1 + exp(x) / nu)
  exp(log_dens)
}

#' Density of Centered log-GED^2 Measurement Noise
#'
#' @param y Numeric vector. Evaluation points (centered: E[eps] = 0).
#' @param nu Numeric. GED shape parameter.
#' @return Density values.
#' @keywords internal
density_eps_ged <- function(y, nu) {
  a <- exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu)))
  mu_bar_g <- (2 / nu) * digamma(1 / nu) + lgamma(1 / nu) - lgamma(3 / nu)
  x <- y + mu_bar_g
  u_abs <- exp(x / 2)
  log_dens <- log(nu) - log(2) - log(a) - lgamma(1 / nu) -
    (u_abs / a)^nu + x / 2
  exp(log_dens)
}


# --------------------------------------------------------------------------- #
# KSC Gaussian Mixture Fitting
# --------------------------------------------------------------------------- #

#' Fit K-Component Gaussian Mixture to Measurement Noise Density
#'
#' Uses EM algorithm to approximate the measurement noise density with a
#' Gaussian mixture. For Gaussian SV, returns the pre-computed KSC (1998)
#' 7-component table.
#'
#' @param distribution Character: \code{"gaussian"}, \code{"student_t"}, or
#'   \code{"ged"}.
#' @param nu Numeric. Shape parameter (ignored for Gaussian).
#' @param K Integer. Number of mixture components. Default 7.
#' @param n_sample Integer. Sample size for EM fitting. Default 500000.
#' @param max_iter Integer. Maximum EM iterations. Default 500.
#' @param tol Numeric. Convergence tolerance. Default 1e-8.
#' @param seed Integer. Random seed. Default 42.
#' @return List with \code{weights}, \code{means}, \code{vars}, \code{KL_div}.
#' @keywords internal
fit_ksc_mixture <- function(distribution = c("gaussian", "student_t", "ged"),
                            nu = NULL, K = 7, n_sample = 100000,
                            max_iter = 500, tol = 1e-8, seed = 42) {
  distribution <- match.arg(distribution)

  # Pre-computed KSC (1998) Table 4 for Gaussian (7 components)
  # These approximate the raw log(chi^2_1) density (mean ≈ -1.27, var ≈ pi^2/2)
  if (distribution == "gaussian" && K == 7) {
    w_ksc <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)
    w_ksc <- w_ksc / sum(w_ksc)  # normalize to sum exactly to 1
    return(list(
      weights = w_ksc,
      means   = c(-11.40039, -5.24321, -9.83726, 1.50746, -0.65098, 0.52478, -2.35859),
      vars    = c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261),
      KL_div  = NA_real_,
      nu      = NULL,
      distribution = "gaussian"
    ))
  }

  # Generate large sample from exact density
  set.seed(seed)
  eps <- switch(distribution,
    gaussian = {
      log(rchisq(n_sample, 1)) - (digamma(0.5) + log(2))
    },
    student_t = {
      if (is.null(nu) || nu <= 0) stop("nu must be positive for Student-t.")
      mu_bar <- digamma(0.5) - digamma(nu / 2) + log(nu)
      log(rf(n_sample, 1, nu)) - mu_bar
    },
    ged = {
      if (is.null(nu) || nu <= 0) stop("nu must be positive for GED.")
      a <- exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu)))
      x_gam <- rgamma(n_sample, shape = 1 / nu, rate = 1)
      u <- sign(runif(n_sample) - 0.5) * a * x_gam^(1 / nu)
      mu_bar_g <- (2 / nu) * digamma(1 / nu) + lgamma(1 / nu) - lgamma(3 / nu)
      log(u^2 + 1e-300) - mu_bar_g
    }
  )

  # EM algorithm for K-component Gaussian mixture
  # Initialize with K-means-style quantile-based initialization
  breaks <- quantile(eps, probs = seq(0, 1, length.out = K + 1))
  m <- numeric(K)
  s2 <- numeric(K)
  q <- rep(1 / K, K)
  for (k in 1:K) {
    idx <- which(eps >= breaks[k] & eps < breaks[k + 1])
    if (k == K) idx <- which(eps >= breaks[k])
    if (length(idx) < 2) idx <- sample(length(eps), max(2, length(eps) %/% K))
    m[k] <- mean(eps[idx])
    s2[k] <- var(eps[idx])
    q[k] <- length(idx) / n_sample
  }

  n <- length(eps)
  for (iter in 1:max_iter) {
    # E-step: compute log responsibilities
    log_tau <- matrix(0, n, K)
    for (k in 1:K) {
      log_tau[, k] <- log(q[k]) + dnorm(eps, mean = m[k], sd = sqrt(s2[k]), log = TRUE)
    }
    # Normalize (log-sum-exp per row)
    max_log <- apply(log_tau, 1, max)
    log_tau_shifted <- log_tau - max_log
    tau <- exp(log_tau_shifted)
    row_sums <- rowSums(tau)
    tau <- tau / row_sums

    # M-step
    q_new <- colMeans(tau)
    m_new <- numeric(K)
    s2_new <- numeric(K)
    for (k in 1:K) {
      Nk <- sum(tau[, k])
      if (Nk < 1e-10) { Nk <- 1e-10 }
      m_new[k] <- sum(tau[, k] * eps) / Nk
      s2_new[k] <- sum(tau[, k] * (eps - m_new[k])^2) / Nk
      if (s2_new[k] < 1e-10) s2_new[k] <- 1e-10
    }

    # Check convergence
    if (max(abs(m_new - m)) < tol && max(abs(q_new - q)) < tol) break
    q <- q_new
    m <- m_new
    s2 <- s2_new
  }

  # Sort by mean
  ord <- order(m)
  q <- q[ord]
  m <- m[ord]
  s2 <- s2[ord]

  list(
    weights = q,
    means   = m,
    vars    = s2,
    KL_div  = NA_real_,
    nu      = nu,
    distribution = distribution
  )
}


# --------------------------------------------------------------------------- #
# Internal: extract filter parameters from model object
# --------------------------------------------------------------------------- #

.get_filter_params <- function(model) {
  phi <- model$phi
  p <- length(phi)
  delta_p <- if (is.null(model$rho) || is.na(model$rho)) 0 else model$rho
  sigma_y <- model$sigy
  sigma_v <- model$sigv

  # Measurement noise variance: distribution-specific sigma_eps^2
  if (inherits(model, "svp_t") && !is.null(model$v) && is.finite(model$v)) {
    sig_eps2 <- psigamma(0.5, 1) + psigamma(model$v / 2, 1)
  } else if (inherits(model, "svp_ged") && !is.null(model$v) && is.finite(model$v)) {
    sig_eps2 <- (2 / model$v)^2 * psigamma(1 / model$v, 1)
  } else {
    sig_eps2 <- (pi^2) / 2
  }

  # Var(z_t) for leverage prediction covariance
  var_zt <- 1.0
  if (inherits(model, "svp_t") && isTRUE(model$leverage) &&
      !is.null(model$v) && is.finite(model$v) && model$v > 2) {
    var_zt <- model$v / (model$v - 2)
  }

  # Measurement intercept mu_bar
  if (inherits(model, "svp_t") && !is.null(model$v) && is.finite(model$v)) {
    mu_bar <- digamma(0.5) - digamma(model$v / 2) + log(model$v)
  } else if (inherits(model, "svp_ged") && !is.null(model$v) && is.finite(model$v)) {
    mu_bar <- (2 / model$v) * digamma(1 / model$v) +
      lgamma(1 / model$v) - lgamma(3 / model$v)
  } else {
    mu_bar <- digamma(0.5) + log(2)
  }
  mu_intercept <- log(sigma_y^2) + mu_bar

  # Distribution name for GMKF
  dist_name <- if (inherits(model, "svp_t")) "student_t"
               else if (inherits(model, "svp_ged")) "ged"
               else "gaussian"

  list(phi = phi, p = p, delta_p = delta_p, sigma_y = sigma_y,
       sigma_v = sigma_v, sig_eps2 = sig_eps2, var_zt = var_zt,
       mu_intercept = mu_intercept, dist_name = dist_name,
       nu = if (!is.null(model$v)) model$v else NULL,
       mu = model$mu)
}


# --------------------------------------------------------------------------- #
# Kalman Filter (CKF) — backward-compatible, enriched output
# --------------------------------------------------------------------------- #

#' Kalman Filter and Smoother for SV(p) Models
#'
#' Applies corrected Kalman filtering (CKF) and Rauch-Tung-Striebel smoothing
#' to estimate the latent log-volatility process from an estimated SV(p) model.
#' Uses distribution-specific measurement noise variance
#' \eqn{\sigma_\varepsilon^2(\nu)}.
#'
#' @param y Numeric vector. Observed returns.
#' @param model An \code{"svp"}, \code{"svp_t"}, or \code{"svp_ged"} object
#'   from \code{\link{svp}}.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#'
#' @return A list with:
#' \describe{
#'   \item{w_estimated}{Filtered log-volatility.}
#'   \item{w_smoothed}{Smoothed log-volatility.}
#'   \item{zt}{Filtered standardized residuals.}
#'   \item{zt_smoothed}{Smoothed standardized residuals.}
#'   \item{P_filtered}{Filtered MSE (first state component).}
#'   \item{P_predicted}{Predicted MSE (first state component).}
#'   \item{xi_filtered}{Full filtered state vectors (p x T matrix).}
#'   \item{xi_smoothed}{Full smoothed state vectors (p x T matrix).}
#'   \item{loglik}{Approximate Gaussian log-likelihood.}
#' }
#'
#' @examples
#' \donttest{
#' sim <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, leverage = TRUE, rho = -0.3)
#' fit <- svp(sim$y, p = 1)
#' kf <- kalman_filter(sim$y, fit)
#' plot(kf$w_smoothed, type = "l", main = "Smoothed log-volatility")
#' }
#'
#' @export
kalman_filter <- function(y, model, del = 1e-10) {
  if (!inherits(model, c("svp", "svp_t", "svp_ged"))) {
    stop("model must be of class 'svp', 'svp_t', or 'svp_ged'.")
  }
  y_vec <- as.numeric(y)
  params <- .get_filter_params(model)

  # Log-squared observations (centered)
  y_star <- log(y_vec^2 + del) - params$mu

  # Build companion matrix F for Lyapunov
  p <- params$p
  F_mat <- matrix(0, p, p)
  F_mat[1, ] <- params$phi
  if (p > 1) for (j in 1:(p - 1)) F_mat[j + 1, j] <- 1

  # Lyapunov initialization: P0 = F P0 F' + sigma_v^2 * r r'
  r_vec <- c(1, rep(0, p - 1))
  Q_init <- params$sigma_v^2 * (r_vec %*% t(r_vec))
  P0 <- tryCatch(
    solve_lyapunov_discrete(F_mat, Q_init),
    error = function(e) diag(p) * params$sigma_v^2  # fallback
  )

  # Call C++ filter
  result <- kalman_filter_cpp(
    y_star = as.numeric(y_star),
    y_raw = y_vec,
    phi = params$phi,
    sigma_y = params$sigma_y,
    sigma_v = params$sigma_v,
    delta_p = params$delta_p,
    sig_eps2 = params$sig_eps2,
    var_zt = params$var_zt,
    P0 = P0
  )

  # Backward compatibility: return w_estimated as matrix
  result$w_estimated <- matrix(result$w_filtered, ncol = 1)
  result$w_smoothed <- matrix(result$w_smoothed, ncol = 1)
  result$zt <- matrix(result$zt, ncol = 1)
  result$zt_smoothed <- matrix(result$zt_smoothed, ncol = 1)
  result
}


# --------------------------------------------------------------------------- #
# filter_svp() — Main user-facing filter function
# --------------------------------------------------------------------------- #

#' Filter Latent Volatility from an Estimated SV(p) Model
#'
#' Applies Kalman filtering (corrected or Gaussian mixture) and RTS smoothing
#' to extract the latent log-volatility process from an estimated SV(p) model.
#'
#' @param object An \code{"svp"}, \code{"svp_t"}, or \code{"svp_ged"} object
#'   from \code{\link{svp}}.
#' @param method Character. Filter method: \code{"corrected"} (default) for
#'   standard Kalman with distribution-specific \eqn{\sigma_\varepsilon^2(\nu)},
#'   \code{"mixture"} for the Gaussian Mixture Kalman Filter (GMKF), or
#'   \code{"particle"} for the Bootstrap Particle Filter (BPF).
#' @param K Integer. Number of mixture components for GMKF. Default 7.
#' @param M Integer. Number of particles for BPF. Default 1000.
#' @param seed Integer. Random seed for BPF. Default 42.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#'
#' @return An object of class \code{"svp_filter"}, a list containing:
#' \describe{
#'   \item{w_filtered}{Filtered log-volatility (T-vector).}
#'   \item{w_smoothed}{Smoothed log-volatility (T-vector).}
#'   \item{zt}{Filtered standardized residuals.}
#'   \item{zt_smoothed}{Smoothed standardized residuals.}
#'   \item{P_filtered}{Filtered MSE of first state component.}
#'   \item{P_predicted}{Predicted MSE of first state component.}
#'   \item{xi_filtered}{Full filtered state vectors (p x T matrix).}
#'   \item{xi_smoothed}{Full smoothed state vectors (p x T matrix).}
#'   \item{loglik}{Approximate log-likelihood.}
#'   \item{method}{Filter method used.}
#'   \item{model}{The input model object.}
#' }
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)
#' fit <- svp(y, p = 1)
#' filt <- filter_svp(fit)
#' plot(filt$w_smoothed, type = "l")
#' }
#'
#' @export
filter_svp <- function(object, method = c("corrected", "mixture", "particle"),
                       K = 7, M = 1000, seed = 42, del = 1e-10) {
  if (!inherits(object, c("svp", "svp_t", "svp_ged"))) {
    stop("object must be of class 'svp', 'svp_t', or 'svp_ged'.")
  }
  method <- match.arg(method)
  y_vec <- as.numeric(object$y)
  params <- .get_filter_params(object)

  # Build companion matrix F
  p <- params$p
  F_mat <- matrix(0, p, p)
  F_mat[1, ] <- params$phi
  if (p > 1) for (j in 1:(p - 1)) F_mat[j + 1, j] <- 1

  # Lyapunov initialization
  r_vec <- c(1, rep(0, p - 1))
  Q_init <- params$sigma_v^2 * (r_vec %*% t(r_vec))
  P0 <- tryCatch(
    solve_lyapunov_discrete(F_mat, Q_init),
    error = function(e) diag(p) * params$sigma_v^2
  )

  if (method == "corrected") {
    # CKF: standard corrected Kalman filter
    y_star <- log(y_vec^2 + del) - params$mu

    result <- kalman_filter_cpp(
      y_star = as.numeric(y_star),
      y_raw = y_vec,
      phi = params$phi,
      sigma_y = params$sigma_y,
      sigma_v = params$sigma_v,
      delta_p = params$delta_p,
      sig_eps2 = params$sig_eps2,
      var_zt = params$var_zt,
      P0 = P0
    )

  } else if (method == "mixture") {
    # GMKF: Gaussian Mixture Kalman Filter
    mixture <- fit_ksc_mixture(params$dist_name, params$nu, K)

    # For GMKF, y_star is raw log(y^2) (NOT centered by mu).
    # The intercept depends on whether mixture means are raw or centered:
    #   - KSC Gaussian table: means ≈ -1.27 (raw log-chi2 density)
    #     → intercept = log(sigma_y^2) only (m_k already encodes E[log(z^2)])
    #   - EM-fitted (Student-t, GED): means ≈ 0 (centered)
    #     → intercept = log(sigma_y^2) + mu_bar (full mu_intercept)
    y_star_raw <- log(y_vec^2 + del)
    mix_mean <- sum(mixture$weights * mixture$means)
    # If mixture mean is far from 0, means are raw → use log(sigy^2) only
    # If near 0, means are centered → use full mu_intercept
    if (abs(mix_mean) > 0.5) {
      mu_intercept_gmkf <- log(params$sigma_y^2)
    } else {
      mu_intercept_gmkf <- params$mu_intercept
    }

    result <- gmkf_filter_cpp(
      y_star = as.numeric(y_star_raw),
      y_raw = y_vec,
      phi = params$phi,
      sigma_y = params$sigma_y,
      sigma_v = params$sigma_v,
      delta_p = params$delta_p,
      var_zt = params$var_zt,
      mix_weights = mixture$weights,
      mix_means = mixture$means,
      mix_vars = mixture$vars,
      mu_intercept = mu_intercept_gmkf,
      P0 = P0
    )
    result$mixture <- mixture

  } else if (method == "particle") {
    # BPF: Bootstrap Particle Filter (C++)
    dist_code <- switch(params$dist_name,
      "gaussian"  = 0L,
      "student_t" = 1L,
      "ged"       = 2L
    )
    nu_val <- if (is.null(params$nu)) 0.0 else params$nu

    result <- particle_filter_svp_cpp(
      y_raw = y_vec,
      phi = params$phi,
      sigma_y = params$sigma_y,
      sigma_v = params$sigma_v,
      nu = nu_val,
      dist_code = dist_code,
      delta = params$delta_p,
      M = as.integer(M),
      seed = as.integer(seed)
    )

    # BPF doesn't produce smoothed states — use filtered as placeholder
    T_obs <- length(y_vec)
    result$w_smoothed <- result$w_filtered  # no smoother for PF
    result$w_predicted <- rep(NA_real_, T_obs)
    result$zt <- y_vec / (params$sigma_y * exp(result$w_filtered / 2))
    result$zt_smoothed <- result$zt
    result$P_predicted <- rep(NA_real_, T_obs)
    # Construct full p x T state matrix from BPF w_filtered
    # Companion state: xi_t = [w_t, w_{t-1}, ..., w_{t-p+1}]'
    xi_mat <- matrix(NA_real_, nrow = p, ncol = T_obs)
    xi_mat[1, ] <- result$w_filtered
    if (p > 1) {
      for (j in 2:p) {
        xi_mat[j, j:T_obs] <- result$w_filtered[1:(T_obs - j + 1)]
      }
    }
    result$xi_filtered <- xi_mat
    result$xi_smoothed <- xi_mat
  }

  # Build output
  p_out <- params$p
  # Full p x p filtered covariance at T (for forecast MSE recursion)
  if (!is.null(result$P_filt_T)) {
    P_filt_T_mat <- matrix(result$P_filt_T, nrow = p_out, ncol = p_out)
  } else {
    # BPF fallback: diagonal approximation
    P_filt_T_mat <- result$P_filtered[length(result$P_filtered)] * diag(p_out)
  }
  out <- list(
    w_filtered  = as.numeric(result$w_filtered),
    w_smoothed  = as.numeric(result$w_smoothed),
    w_predicted = as.numeric(result$w_predicted),
    zt          = as.numeric(result$zt),
    zt_smoothed = as.numeric(result$zt_smoothed),
    P_filtered  = as.numeric(result$P_filtered),
    P_predicted = as.numeric(result$P_predicted),
    P_filt_T    = P_filt_T_mat,
    xi_filtered = result$xi_filtered,
    xi_smoothed = result$xi_smoothed,
    loglik      = result$loglik,
    method      = method,
    model       = object
  )
  if (!is.null(result$mixture)) out$mixture <- result$mixture
  if (!is.null(result$ESS)) out$ESS <- as.numeric(result$ESS)
  class(out) <- "svp_filter"
  out
}
