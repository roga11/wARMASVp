# =========================================================================== #
# Internal estimation helpers called by svp() in estim.R
# =========================================================================== #


# =========================================================================== #
# Correction factors and helper functions for leverage under heavy tails
# =========================================================================== #

#' Correction factor C_t(nu) for Student-t leverage estimation
#'
#' Under scale mixture \eqn{u_t = z_t \lambda_t^{-1/2}},
#' \eqn{C_t(\nu) = [E[\lambda^{-1/2}]]^2 = (\nu/2) [\Gamma((\nu-1)/2) / \Gamma(\nu/2)]^2}.
#' Exact, parameter-free.
#' @param nu Degrees of freedom (nu > 1)
#' @return C_t(nu), always > 1 for finite nu (approaches 1 as nu -> Inf)
#' @keywords internal
correction_factor_t <- function(nu) {
  exp(log(nu / 2) + 2 * (lgamma((nu - 1) / 2) - lgamma(nu / 2)))
}

#' E[|u|] for standardized GED(nu) with Var = 1
#' Closed form: sqrt(Gamma(1/nu)/Gamma(3/nu)) * Gamma(2/nu) / Gamma(1/nu)
#' @param nu GED shape parameter (nu > 0)
#' @return E[|u|]
#' @keywords internal
ged_E_abs_u <- function(nu) {
  exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu)) + lgamma(2 / nu) - lgamma(1 / nu))
}

#' Quantile function for standardized GED(nu) with Var = 1
#' Uses the relationship between GED CDF and incomplete gamma function.
#' @param p Probability (0 < p < 1). Clamped to [1e-15, 1-1e-15] internally.
#' @param nu GED shape parameter (nu > 0)
#' @return Quantile value
#' @keywords internal
qged_std <- function(p, nu) {
  p <- pmax(1e-15, pmin(1 - 1e-15, p))
  a <- exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu)))
  # Handle p = 0.5 separately (returns 0) to avoid qgamma(0, ...) NaN warnings
  result <- numeric(length(p))
  hi <- p > 0.5
  lo <- p < 0.5
  if (any(hi)) {
    result[hi] <- a * qgamma(2 * (p[hi] - 0.5), shape = 1 / nu)^(1 / nu)
  }
  if (any(lo)) {
    result[lo] <- -a * qgamma(2 * (0.5 - p[lo]), shape = 1 / nu)^(1 / nu)
  }
  result
}

#' Gauss-Hermite quadrature nodes and weights for N(0,1) integration
#' Computes nodes z_i and weights w_i such that sum(w_i * f(z_i)) approximates E[f(Z)]
#' where Z ~ N(0,1). Uses the Golub-Welsch algorithm.
#' @param n Number of quadrature points (default 200)
#' @return List with components nodes and weights
#' @keywords internal
.gauss_hermite_normal <- function(n = 200L) {
  # Golub-Welsch: eigenvalues of tridiagonal Jacobi matrix for Hermite polynomials
  i <- seq_len(n)
  # Recurrence coefficients for (physicist's) Hermite: a_i = 0, b_i = sqrt(i/2)
  b <- sqrt(i[-n] / 2)
  # Tridiagonal matrix
  J <- matrix(0, n, n)
  diag(J) <- 0
  for (k in seq_len(n - 1L)) {
    J[k, k + 1L] <- b[k]
    J[k + 1L, k] <- b[k]
  }
  eig <- eigen(J, symmetric = TRUE)
  # Physicist's Hermite nodes: multiply by sqrt(2) for N(0,1)
  nodes <- eig$values * sqrt(2)
  # Weights: first component squared, normalized for N(0,1)
  weights <- eig$vectors[1, ]^2
  # Sort by nodes
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
}

# Package-level cached GH quadrature (lazy initialization)
.gh_cache <- new.env(parent = emptyenv())

#' Get cached Gauss-Hermite nodes/weights for N(0,1)
#' @param n Number of nodes (default 200)
#' @return List with nodes and weights
#' @keywords internal
.get_gh <- function(n = 200L) {
  key <- paste0("gh_", n)
  if (is.null(.gh_cache[[key]])) {
    .gh_cache[[key]] <- .gauss_hermite_normal(n)
  }
  .gh_cache[[key]]
}

#' Approximate correction factor C_g(nu) for GED leverage (diagnostic use only)
#'
#' \eqn{C_g(\nu) = E[|u|] \cdot \textrm{Cov}(z, F_{GED}^{-1}(\Phi(z))) / \sqrt{2/\pi}}.
#' This is a first-order approximation; estimation uses the exact implicit equation.
#' @param nu GED shape parameter (nu > 0)
#' @param n_nodes Number of GH quadrature nodes (default 200)
#' @return Approximate C_g(nu)
#' @keywords internal
correction_factor_ged_approx <- function(nu, n_nodes = 200L) {
  E_abs_u <- ged_E_abs_u(nu)
  gh <- .get_gh(n_nodes)
  u_vals <- qged_std(stats::pnorm(gh$nodes), nu)
  cov_zu <- sum(gh$weights * gh$nodes * u_vals)
  E_abs_u * cov_zu / sqrt(2 / pi)
}


# =========================================================================== #
# Unified leverage estimation — works for all distributions
# =========================================================================== #

#' Compute EH cross-moment for leverage estimation
#'
#' EH = demeaned sample cross-moment \eqn{(1/(T-2)) \sum(|y_t| - \bar{|y|})(y_{t-1} - \bar{y})}.
#' @param y Numeric vector of observations
#' @param rho_type "pearson" or "kendall"
#' @return EH value
#' @keywords internal
.compute_EH <- function(y, rho_type = "pearson") {
  N <- length(y)
  yabs <- abs(y)
  muu <- mean(y[1:(N - 1)])
  mua <- mean(yabs[2:N])
  if (rho_type == "kendall") {
    tau <- kendall_corr(yabs[2:N] - mua, y[1:(N - 1)] - muu)
    EH <- tau * sqrt(stats::var(yabs[2:N]) * stats::var(y[1:(N - 1)]))
  } else {
    EH <- sum((yabs[2:N] - mua) * (y[1:(N - 1)] - muu)) / (N - 2)
  }
  as.numeric(EH)
}

#' Add leverage estimation to a fitted SV(p) model (post-processing)
#' Works for all error distributions: Gaussian, Student-t, GED.
#' - Gaussian: closed form (C_F = 1)
#' - Student-t: closed form with C_t(nu) correction
#' - GED: exact 1D root-finding via uniroot + Gauss-Hermite quadrature
#' @param out Model object from .svp_gaussian/.svp_t/.svp_ged (without leverage)
#' @param y Numeric vector of observations
#' @param p AR order
#' @param rho_type "pearson", "kendall", or "both"
#' @param del Small constant for log transform
#' @param trunc_lev Logical: truncate rho to [-0.999, 0.999]
#' @param wDecay Logical: decaying weights (passed from original estimation)
#' @param errorType "Gaussian", "Student-t", or "GED"
#' @return Updated model object with leverage fields added
#' @keywords internal
.add_leverage <- function(out, y, p, rho_type, del, trunc_lev, wDecay, errorType) {
  y <- as.numeric(y)
  phi <- as.numeric(out$phi)
  sigy <- out$sigy
  sigv <- out$sigv
  nu <- if (!is.null(out$v)) out$v else NULL

  # --- Compute gammatilde (general-p, distribution-free) ---
  rho_w <- as.numeric(stats::ARMAacf(ar = phi, lag.max = p)[-1])
  gamma_w0 <- sigv^2 / (1 - sum(phi * rho_w))
  gammatilde <- gamma_w0 * (1 + rho_w[1])

  # --- Guard against degenerate cases ---
  if (sigv < 1e-10 || sigy < 1e-10) {
    warning("Leverage estimation failed: sigma_v or sigma_y near zero.")
    out$rho <- NA_real_
    out$gammatilde <- gammatilde
    out$leverage <- TRUE
    out$rho_type <- rho_type
    out$trunc_lev <- trunc_lev
    out$theta <- c(out$theta, NA_real_)
    return(out)
  }

  # --- Helper to compute delta from EH ---
  .compute_delta <- function(EH_val) {
    if (errorType == "Gaussian") {
      # Closed form: delta = sqrt(2pi) * EH / (sigv * sigy^2) * exp(-gammatilde/4)
      delta <- (sqrt(2 * pi) * EH_val) / (sigv * sigy^2) * exp(-0.25 * gammatilde)

    } else if (errorType == "Student-t") {
      # Closed form with C_t(nu) correction
      Ct <- correction_factor_t(nu)
      delta <- (sqrt(2 * pi) * EH_val) / (sigv * sigy^2 * Ct) * exp(-0.25 * gammatilde)

    } else if (errorType == "GED") {
      # Exact 1D root-finding via C++:
      # Solve E[g(z + delta*sigv/2)] = EH / (sigy^2 * E[|u|] * exp(gammatilde/4))
      E_abs_u <- ged_E_abs_u(nu)
      target <- EH_val / (sigy^2 * E_abs_u * exp(gammatilde / 4))
      gh <- .get_gh(200L)
      delta <- find_delta_ged_cpp(target, sigv, nu, gh$nodes, gh$weights)
      if (is.na(delta)) {
        warning("GED leverage root-finding failed.")
      }
    }
    delta
  }

  # --- Compute leverage for requested rho_type(s) ---
  if (rho_type == "pearson") {
    EH <- .compute_EH(y, "pearson")
    rho_val <- .compute_delta(EH)
    if (trunc_lev && !is.na(rho_val)) rho_val <- max(-0.999, min(0.999, rho_val))
  } else if (rho_type == "kendall") {
    EH <- .compute_EH(y, "kendall")
    rho_val <- .compute_delta(EH)
    if (trunc_lev && !is.na(rho_val)) rho_val <- max(-0.999, min(0.999, rho_val))
  } else if (rho_type == "both") {
    EH_pea <- .compute_EH(y, "pearson")
    EH_ken <- .compute_EH(y, "kendall")
    rho_val <- .compute_delta(EH_pea)
    rho_ken <- .compute_delta(EH_ken)
    if (trunc_lev && !is.na(rho_val)) rho_val <- max(-0.999, min(0.999, rho_val))
    if (trunc_lev && !is.na(rho_ken)) rho_ken <- max(-0.999, min(0.999, rho_ken))
    out$rho_kendall <- rho_ken
    out$rho_pearson <- rho_val
  } else {
    rho_val <- NA_real_
  }

  # --- Update model object ---
  out$rho <- rho_val
  out$gammatilde <- gammatilde
  out$leverage <- TRUE
  out$rho_type <- rho_type
  out$trunc_lev <- trunc_lev
  out$theta <- c(out$theta, rho_val)

  # Store correction factor for diagnostics
  if (errorType == "Student-t") {
    out$CF <- correction_factor_t(nu)
  } else if (errorType == "GED") {
    out$CF <- correction_factor_ged_approx(nu)  # approximate, for diagnostic reporting
  } else {
    out$CF <- 1
  }

  return(out)
}


# --- Gaussian SV(p) estimation (with optional leverage) ---
.svp_gaussian <- function(y, p, J, leverage, rho_type, del, trunc_lev, wDecay,
                          sigvMethod = "factored") {
  y <- as.numeric(y)
  if (length(y) < 2 * p + J) {
    stop("Time series too short for the given p and J.")
  }
  if (!rho_type %in% c("pearson", "kendall", "both", "none")) {
    stop("rho_type must be one of 'pearson', 'kendall', 'both', 'none'.")
  }
  # Step 1: phi estimation (C++, distribution-free)
  para <- svpCpp_nolev(y, p, J, del, wDecay)
  # Step 2: sigv override by method
  # C++ computes hybrid (factored for p=1, direct for p>=2) by default
  ly2 <- log(y^2 + del)
  ys_g <- ly2 - mean(ly2)
  N_g <- length(ly2)
  var_ly2 <- stats::var(as.numeric(ly2))
  sig_eps2 <- (pi^2) / 2
  phi_g <- as.numeric(para$phi)
  if (sigvMethod == "direct") {
    gam_vec_g <- numeric(p)
    for (k in seq_len(p)) {
      gam_vec_g[k] <- (1 / (N_g - 1)) * sum(ys_g[(k + 1):N_g] * ys_g[1:(N_g - k)])
    }
    sv2 <- var_ly2 - sum(phi_g * gam_vec_g) - sig_eps2
    para$sigv <- sqrt(abs(sv2))
  } else if (sigvMethod == "factored") {
    rho_w_p <- as.numeric(stats::ARMAacf(ar = phi_g, lag.max = p)[-1])
    sv2 <- (var_ly2 - sig_eps2) * (1 - sum(phi_g * rho_w_p))
    para$sigv <- sqrt(abs(sv2))
  }
  # else "hybrid": keep C++ result as-is
  para$theta <- c(para$phi, para$sigy, para$sigv)
  # Metadata
  para$rho <- NA_real_
  para$gammatilde <- NA_real_
  para$y <- y
  para$p <- p
  para$J <- J
  para$leverage <- FALSE
  para$del <- del
  para$wDecay <- wDecay
  para$trunc_lev <- trunc_lev
  para$rho_type <- NA_character_
  para$sigvMethod <- sigvMethod
  class(para) <- "svp"
  # Step 3: leverage post-processing (if requested)
  if (leverage) {
    para <- .add_leverage(para, y, p, rho_type, del, trunc_lev, wDecay, "Gaussian")
  }
  return(para)
}

# --- Student-t SV(p) estimation ---
.svp_t <- function(y, p, J, del, wDecay, logNu,
                   sigvMethod = "factored", winsorize_eps = 0) {
  y <- as.matrix(as.numeric(y))
  N <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  # Estimate phi via W-ARMA-SV (distribution-free)
  para_base <- svpCpp_nolev(as.numeric(y), p, J, del, wDecay)
  phi_reg <- as.numeric(para_base$phi)
  var_ly2 <- stats::var(as.numeric(ly2))
  # Estimate sigma_eps^2
  # winsorize_eps: 0 or FALSE = off; TRUE = use J; integer > 0 = use that as J_w
  J_w <- 0L
  if (is.logical(winsorize_eps)) {
    if (winsorize_eps) J_w <- J
  } else if (is.numeric(winsorize_eps) && winsorize_eps > 0) {
    J_w <- as.integer(winsorize_eps)
  }
  if (J_w > 0L) {
    # OLS winsorized estimator (SVHT Eq. 69): average gamma_w(0) across J_w lags
    rho_w_all <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = J_w)[-1])
    gam_all <- numeric(J_w)
    for (k in seq_len(J_w)) {
      gam_all[k] <- (1 / (N - 1)) * as.numeric(
        t(ys[(k + 1):N, , drop = FALSE]) %*% ys[1:(N - k), , drop = FALSE])
    }
    gamma_w0_ols <- sum(gam_all * rho_w_all) / sum(rho_w_all^2)
    se2 <- var_ly2 - gamma_w0_ols
  } else {
    # Single-lag estimator (SVHT Eq. 18): sigma_eps^2 = gamma(0) - gamma(1)/rho_w(1)
    rho_w1 <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = 1)[2])
    gam1 <- (1 / (N - 1)) * as.numeric(
      t(ys[2:N, , drop = FALSE]) %*% ys[1:(N - 1), , drop = FALSE])
    se2 <- var_ly2 - gam1 / rho_w1
  }
  se2b <- as.numeric(se2) - psigamma(0.5, 1)
  # Estimate nu via root-finding on trigamma: psigamma(nu/2, 1) = se2b (C++)
  nu_lower <- 2.01
  nu_upper <- 500
  nuh <- find_nu_t_cpp(se2b, nu_lower, nu_upper, logNu)
  if (nuh >= nu_upper) {
    warning("Estimated nu at upper boundary (", nu_upper,
            "); tails indistinguishable from Gaussian.")
  } else if (nuh <= nu_lower) {
    warning("Estimated nu at lower boundary (", nu_lower,
            "); extremely heavy tails.")
  }
  # Estimate sigv
  var_log_sq_t <- psigamma(0.5, 1) + psigamma(nuh / 2, 1)
  # Compute sample autocovariances at lags 1,...,p (needed by direct/hybrid for p>=2)
  gam_vec <- numeric(p)
  for (k in seq_len(p)) {
    gam_vec[k] <- (1 / (N - 1)) * as.numeric(
      t(ys[(k + 1):N, , drop = FALSE]) %*% ys[1:(N - k), , drop = FALSE])
  }
  if (sigvMethod == "direct") {
    sv_reg <- sqrt(abs(var_ly2 - sum(phi_reg * gam_vec) - var_log_sq_t))
  } else if (sigvMethod == "factored") {
    rho_w_p <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = p)[-1])
    sv_reg <- sqrt(abs((var_ly2 - var_log_sq_t) * (1 - sum(phi_reg * rho_w_p))))
  } else {
    # "hybrid": AD2021 for p=1, direct for p>=2
    if (p == 1L) {
      sv_reg <- sqrt(abs((1 - phi_reg^2) * (var_ly2 - var_log_sq_t)))
    } else {
      sv_reg <- sqrt(abs(var_ly2 - sum(phi_reg * gam_vec) - var_log_sq_t))
    }
  }
  # Estimate sigy
  mu_log_sq_t <- psigamma(0.5, 0) - psigamma(nuh / 2, 0) + log(nuh)
  if (p == 1L) {
    sigy <- sqrt(exp(mean(log(y^2)) - mu_log_sq_t))
  } else {
    sigy <- sqrt(exp(mean(log(y^2 + del)) - mu_log_sq_t))
  }
  theta <- c(phi_reg, sigy, sv_reg, nuh)
  out <- list(mu = mu, phi = phi_reg, sigv = as.numeric(sv_reg),
              sigy = as.numeric(sigy), v = nuh, theta = theta,
              y = as.numeric(y), J = J, p = as.integer(p), del = del,
              wDecay = wDecay, logNu = logNu,
              sigvMethod = sigvMethod, winsorize_eps = winsorize_eps,
              rho = NA_real_, gammatilde = NA_real_,
              leverage = FALSE, rho_type = NA_character_,
              trunc_lev = TRUE,
              nonstationary_ind = isTRUE(para_base$nonstationary_ind))
  class(out) <- "svp_t"
  return(out)
}

# --- GED SV(p) estimation ---
.svp_ged <- function(y, p, J, del, wDecay,
                     sigvMethod = "factored", winsorize_eps = 0) {
  y <- as.matrix(as.numeric(y))
  N <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  # Estimate phi (distribution-free)
  para_base <- svpCpp_nolev(as.numeric(y), p, J, del, wDecay)
  phi_reg <- as.numeric(para_base$phi)
  var_ly2 <- stats::var(as.numeric(ly2))
  # Estimate sigma_eps^2
  # winsorize_eps: 0 or FALSE = off; TRUE = use J; integer > 0 = use that as J_w
  J_w <- 0L
  if (is.logical(winsorize_eps)) {
    if (winsorize_eps) J_w <- J
  } else if (is.numeric(winsorize_eps) && winsorize_eps > 0) {
    J_w <- as.integer(winsorize_eps)
  }
  if (J_w > 0L) {
    # OLS winsorized estimator (SVHT Eq. 69)
    rho_w_all <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = J_w)[-1])
    gam_all <- numeric(J_w)
    for (k in seq_len(J_w)) {
      gam_all[k] <- (1 / (N - 1)) * as.numeric(
        t(ys[(k + 1):N, , drop = FALSE]) %*% ys[1:(N - k), , drop = FALSE])
    }
    gamma_w0_ols <- sum(gam_all * rho_w_all) / sum(rho_w_all^2)
    se2 <- var_ly2 - gamma_w0_ols
  } else {
    # Single-lag estimator (SVHT Eq. 18)
    rho_w1 <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = 1)[2])
    gam1 <- (1 / (N - 1)) * as.numeric(
      t(ys[2:N, , drop = FALSE]) %*% ys[1:(N - 1), , drop = FALSE])
    se2 <- as.numeric(var_ly2 - gam1 / rho_w1)
  }
  # Estimate nu: solve (2/nu)^2 * psigamma(1/nu, 1) = se2 (C++)
  m1 <- sum(abs(ly2)) / N
  m2 <- sqrt(sum((abs(ly2) - m1)^2) / N)
  x0 <- m1 / m2
  lower <- max(1e-6, x0 / 5)
  upper <- x0 * 5
  ged_nuh <- find_nu_ged_cpp(as.numeric(se2), lower, upper)
  if (ged_nuh >= 20) {
    warning("Estimated GED nu at upper boundary (20); tails indistinguishable from Gaussian.")
  } else if (ged_nuh <= 0.1) {
    warning("Estimated GED nu at lower boundary (0.1); extremely heavy tails.")
  }
  # Estimate sigv
  var_log_sq_ged <- ((2 / ged_nuh)^2) * psigamma(1 / ged_nuh, 1)
  gam_vec <- numeric(p)
  for (k in seq_len(p)) {
    gam_vec[k] <- (1 / (N - 1)) * as.numeric(
      t(ys[(k + 1):N, , drop = FALSE]) %*% ys[1:(N - k), , drop = FALSE])
  }
  if (sigvMethod == "direct") {
    sv_reg <- sqrt(abs(var_ly2 - sum(phi_reg * gam_vec) - var_log_sq_ged))
  } else if (sigvMethod == "factored") {
    rho_w_p <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = p)[-1])
    sv_reg <- sqrt(abs((var_ly2 - var_log_sq_ged) * (1 - sum(phi_reg * rho_w_p))))
  } else {
    if (p == 1L) {
      sv_reg <- sqrt(abs((1 - phi_reg^2) * (var_ly2 - var_log_sq_ged)))
    } else {
      sv_reg <- sqrt(abs(var_ly2 - sum(phi_reg * gam_vec) - var_log_sq_ged))
    }
  }
  # Estimate sigy
  mu_log_sq_ged <- (2 / ged_nuh) * psigamma(1 / ged_nuh, 0) +
    log(gamma(1 / ged_nuh)) - log(gamma(3 / ged_nuh))
  if (p == 1L) {
    sigy <- sqrt(exp(mean(log(y^2)) - mu_log_sq_ged))
  } else {
    sigy <- sqrt(exp(mean(log(y^2 + del)) - mu_log_sq_ged))
  }
  theta <- c(phi_reg, sigy, sv_reg, ged_nuh)
  out <- list(mu = mu, phi = phi_reg, sigv = as.numeric(sv_reg),
              sigy = as.numeric(sigy), v = ged_nuh, theta = theta,
              y = as.numeric(y), J = J, p = as.integer(p), del = del,
              wDecay = wDecay,
              sigvMethod = sigvMethod, winsorize_eps = winsorize_eps,
              rho = NA_real_, gammatilde = NA_real_,
              leverage = FALSE, rho_type = NA_character_,
              trunc_lev = TRUE,
              nonstationary_ind = isTRUE(para_base$nonstationary_ind))
  class(out) <- "svp_ged"
  return(out)
}








# --- Internal SE helpers ---

.svpSE_gaussian <- function(object, n_sim, alpha, burnin) {
  p <- object$p
  Tsize <- length(object$y)
  has_lev <- isTRUE(object$leverage) && !is.na(object$rho)
  wDecay <- if (is.null(object$wDecay)) FALSE else object$wDecay
  trunc_lev <- if (is.null(object$trunc_lev)) TRUE else object$trunc_lev
  n_params <- if (has_lev) p + 3 else p + 2
  betamat <- matrix(0, n_sim, n_params)
  if (has_lev) {
    betasim <- c(object$phi, object$sigy, object$sigv, object$rho)
    rho_type <- if (is.null(object$rho_type)) "pearson" else object$rho_type
    if (rho_type == "both") rho_type <- "pearson"
  } else {
    betasim <- c(object$phi, object$sigy, object$sigv)
    rho_type <- "pearson"
  }
  xn <- 1
  while (xn <= n_sim) {
    if (has_lev) {
      u_out <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                       sigv = object$sigv, leverage = TRUE,
                       rho = object$rho, burnin = burnin)
      u <- u_out$y
    } else {
      u <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                   sigv = object$sigv, burnin = burnin)$y
    }
    sigvMethod_g <- if (is.null(object$sigvMethod)) "factored" else object$sigvMethod
    out_tmp <- tryCatch(
      svp(u, p = p, J = object$J, leverage = has_lev,
          rho_type = rho_type, del = object$del, trunc_lev = trunc_lev,
          wDecay = wDecay, sigvMethod = sigvMethod_g),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && !isTRUE(out_tmp$nonstationary_ind)) {
      if (has_lev) {
        betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv, out_tmp$rho)
      } else {
        betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv)
      }
      xn <- xn + 1
    }
  }
  .compute_se_ci(betamat, betasim, n_sim, n_params, alpha)
}

.svpSE_t <- function(object, n_sim, alpha, burnin, logNu) {
  Tsize <- length(object$y)
  p <- object$p
  has_lev <- isTRUE(object$leverage) && !is.na(object$rho)
  wDecay <- if (is.null(object$wDecay)) FALSE else object$wDecay
  trunc_lev <- if (is.null(object$trunc_lev)) TRUE else object$trunc_lev
  sigvMethod_t <- if (is.null(object$sigvMethod)) "factored" else object$sigvMethod
  winsorize_eps_t <- if (is.null(object$winsorize_eps)) FALSE else object$winsorize_eps
  rho_type <- if (is.null(object$rho_type) || is.na(object$rho_type)) "pearson" else object$rho_type
  if (rho_type == "both") rho_type <- "pearson"
  n_params <- if (has_lev) length(object$phi) + 4 else length(object$phi) + 3
  betamat <- matrix(0, n_sim, n_params)
  if (has_lev) {
    betasim <- c(object$phi, object$sigy, object$sigv, object$v, object$rho)
  } else {
    betasim <- c(object$phi, object$sigy, object$sigv, object$v)
  }
  xn <- 1
  while (xn <= n_sim) {
    if (has_lev) {
      u_out <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                       sigv = object$sigv, errorType = "Student-t",
                       leverage = TRUE, rho = object$rho,
                       nu = object$v, burnin = burnin)
      u <- u_out$y
    } else {
      u <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                   sigv = object$sigv, errorType = "Student-t",
                   nu = object$v, burnin = burnin)$y
    }
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = object$J,
          errorType = "Student-t", leverage = has_lev,
          rho_type = rho_type, del = object$del, logNu = logNu,
          trunc_lev = trunc_lev, wDecay = wDecay,
          sigvMethod = sigvMethod_t, winsorize_eps = winsorize_eps_t),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && is.finite(out_tmp$v)) {
      if (has_lev && (is.na(out_tmp$rho) || !is.finite(out_tmp$rho))) next
      if (has_lev) {
        betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv,
                           out_tmp$v, out_tmp$rho)
      } else {
        betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv, out_tmp$v)
      }
      xn <- xn + 1
    }
  }
  .compute_se_ci(betamat, betasim, n_sim, n_params, alpha)
}

.svpSE_ged <- function(object, n_sim, alpha, burnin) {
  Tsize <- length(object$y)
  p <- object$p
  has_lev <- isTRUE(object$leverage) && !is.na(object$rho)
  wDecay <- if (is.null(object$wDecay)) FALSE else object$wDecay
  trunc_lev <- if (is.null(object$trunc_lev)) TRUE else object$trunc_lev
  sigvMethod_ged <- if (is.null(object$sigvMethod)) "factored" else object$sigvMethod
  winsorize_eps_ged <- if (is.null(object$winsorize_eps)) FALSE else object$winsorize_eps
  rho_type <- if (is.null(object$rho_type) || is.na(object$rho_type)) "pearson" else object$rho_type
  if (rho_type == "both") rho_type <- "pearson"
  n_params <- if (has_lev) length(object$phi) + 4 else length(object$phi) + 3
  betamat <- matrix(0, n_sim, n_params)
  if (has_lev) {
    betasim <- c(object$phi, object$sigy, object$sigv, object$v, object$rho)
  } else {
    betasim <- c(object$phi, object$sigy, object$sigv, object$v)
  }
  xn <- 1
  while (xn <= n_sim) {
    if (has_lev) {
      u_out <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                       sigv = object$sigv, errorType = "GED",
                       leverage = TRUE, rho = object$rho,
                       nu = object$v, burnin = burnin)
      u <- u_out$y
    } else {
      u <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                   sigv = object$sigv, errorType = "GED",
                   nu = object$v, burnin = burnin)$y
    }
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = object$J,
          errorType = "GED", leverage = has_lev,
          rho_type = rho_type, del = object$del,
          trunc_lev = trunc_lev, wDecay = wDecay,
          sigvMethod = sigvMethod_ged, winsorize_eps = winsorize_eps_ged),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && is.finite(out_tmp$v)) {
      if (has_lev && (is.na(out_tmp$rho) || !is.finite(out_tmp$rho))) next
      if (has_lev) {
        betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv,
                           out_tmp$v, out_tmp$rho)
      } else {
        betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv, out_tmp$v)
      }
      xn <- xn + 1
    }
  }
  .compute_se_ci(betamat, betasim, n_sim, n_params, alpha)
}

# Internal helper for SE/CI computation
#' @keywords internal
.compute_se_ci <- function(betamat, betasim, n_sim, n_params, alpha) {
  betasim_mat <- matrix(betasim, nrow = n_sim, ncol = n_params, byrow = TRUE)
  SEsim0 <- sqrt(colSums((betamat - betasim_mat)^2) / (n_sim - n_params))
  SEsim <- sqrt(colSums((betamat - matrix(colMeans(betamat), n_sim, n_params, byrow = TRUE))^2) / (n_sim - n_params))
  ISEconservative <- numeric(n_params)
  ISEliberal <- numeric(n_params)
  CI <- matrix(0, n_params, 2)
  for (xp in seq_len(n_params)) {
    sorted <- sort(betamat[, xp])
    cl <- sorted[max(1, round(n_sim * (alpha / 2)))]
    ch <- sorted[round(n_sim * (1 - alpha / 2))]
    z_alpha <- stats::qnorm(1 - alpha / 2)
    sig_L <- (betasim[xp] - cl) / z_alpha
    sig_H <- (ch - betasim[xp]) / z_alpha
    ISEconservative[xp] <- min(abs(sig_L), abs(sig_H))
    ISEliberal[xp] <- max(abs(sig_L), abs(sig_H))
    CI[xp, 1] <- cl
    CI[xp, 2] <- ch
  }
  list(CI = t(CI), SEsim0 = SEsim0, SEsim = SEsim,
       ISEconservative = ISEconservative, ISEliberal = ISEliberal,
       thetamat = betamat)
}


# =========================================================================== #
# Internal helpers for hypothesis testing
# =========================================================================== #

# --- MMC p-value functions (return negative for minimization) ---

# MMC p-value for AR test
.mmc_pval_ar <- function(theta, y, p_null, p_alt, j, N, s0, ini,
                         del = 1e-10, wDecay = FALSE, Bartlett = FALSE,
                         sigvMethod = "factored",
                         errorType = "Gaussian",
                         nu_fixed = NA_real_,
                         logNu = TRUE,
                         winsorize_eps = 0,
                         innovations = NULL) {
  Tsize <- length(y)
  phi_null <- theta[1:p_null]
  sigy_null <- theta[p_null + 1]
  sigv_null <- theta[p_null + 2]
  stationary <- all(Mod(polyroot(c(1, -phi_null))) > 1)
  if (!stationary || sigy_null <= 0 || sigv_null <= 0) {
    return(9999999999999)
  }
  betasim_null <- c(phi_null, sigy_null, sigv_null)
  sN <- .simnull_ar(betasim_null, p_null, p_alt, j, Tsize, N, ini,
                    del, wDecay, Bartlett, sigvMethod = sigvMethod,
                    errorType = errorType, nu_fixed = nu_fixed,
                    logNu = logNu, winsorize_eps = winsorize_eps,
                    innovations = innovations)
  pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  return(pval)
}

# MMC p-value for leverage test (identity Amat)
# S0 is kept fixed per Dufour (2006, eq 4.22)
.mmc_pval_lev <- function(theta, y, j, N, s0, mdl_alt, rho_null, ini,
                          Amat, rho_type, del = 1e-10, Bartlett = FALSE,
                          wDecay = FALSE, trunc_lev = TRUE,
                          logNu = FALSE, sigvMethod = "factored",
                          winsorize_eps = 0,
                          innovations = NULL) {
  p <- length(mdl_alt$phi)
  Tsize <- length(as.numeric(y))
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (!stationary || theta[p + 1] <= 0 || theta[p + 2] <= 0) {
    return(9999999999999)
  }
  betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], rho_null)
  sN <- .simnull(betasim_null, rho_null, p, j, Tsize, N, ini,
                 Amat, rho_type, del, wDecay = wDecay,
                 trunc_lev = trunc_lev, sigvMethod = sigvMethod,
                 innovations = innovations)
  pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  return(pval)
}

# MMC p-value for leverage test (Bartlett kernel Amat)
# S0 is kept fixed per Dufour (2006, eq 4.22)
.mmc_pval_lev_Amat <- function(theta, y, j, N, s0, mdl_alt, rho_null, ini, innovations = NULL,
                               Amat = NULL, rho_type, del = 1e-10, Bartlett = TRUE,
                               wDecay = FALSE, trunc_lev = TRUE,
                               logNu = FALSE, sigvMethod = "factored",
                               winsorize_eps = 0) {
  p <- length(mdl_alt$phi)
  Tsize <- length(as.numeric(y))
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (!stationary || theta[p + 1] <= 0 || theta[p + 2] <= 0) {
    return(9999999999999)
  }
  betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], rho_null)
  sN <- .simnull_Amat(betasim_null, rho_null, p, j, Tsize, N, ini,
                      rho_type, del, Bartlett, wDecay = wDecay,
                      trunc_lev = trunc_lev, sigvMethod = sigvMethod,
                      innovations = innovations)
  pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  return(pval)
}

# MMC p-value for Student-t leverage test (returns negative for minimization)
# theta = (phi_1,...,phi_p, sigy, sigv, nu) — nuisance under H0: rho = rho_null
# S0 is kept fixed per Dufour (2006, eq 4.22)
.mmc_pval_lev_t <- function(theta, y, j, N, s0, mdl_alt, rho_null, ini, innovations = NULL,
                             Amat, rho_type, del = 1e-10, Bartlett = FALSE,
                             wDecay = FALSE, trunc_lev = TRUE,
                             logNu = FALSE, sigvMethod = "factored",
                             winsorize_eps = FALSE) {
  p <- length(mdl_alt$phi)
  Tsize <- length(as.numeric(y))
  phi_cand <- theta[1:p]
  sigy_cand <- theta[p + 1]
  sigv_cand <- theta[p + 2]
  nu_cand <- theta[p + 3]
  stationary <- all(Mod(polyroot(c(1, -phi_cand))) > 1)
  if (!stationary || sigy_cand <= 0 || sigv_cand <= 0 || nu_cand <= 2) {
    return(9999999999999)
  }
  betasim_null <- c(phi_cand, sigy_cand, sigv_cand, nu_cand, rho_null)
  if (isTRUE(Bartlett)) {
    sN <- .simnull_lev_t_Amat(betasim_null, rho_null, p, j, Tsize, N, ini,
                               rho_type, del, TRUE, wDecay = wDecay,
                               trunc_lev = trunc_lev, logNu = logNu,
                               sigvMethod = sigvMethod,
                               winsorize_eps = winsorize_eps,
                               innovations = innovations)
  } else {
    sN <- .simnull_lev_t(betasim_null, rho_null, p, j, Tsize, N, ini,
                          Amat, rho_type, del, wDecay = wDecay,
                          trunc_lev = trunc_lev, logNu = logNu,
                          sigvMethod = sigvMethod,
                          winsorize_eps = winsorize_eps,
                          innovations = innovations)
  }
  pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  return(pval)
}

# MMC p-value for GED leverage test (returns negative for minimization)
# theta = (phi_1,...,phi_p, sigy, sigv, nu) — nuisance under H0: rho = rho_null
# S0 is kept fixed per Dufour (2006, eq 4.22)
.mmc_pval_lev_ged <- function(theta, y, j, N, s0, mdl_alt, rho_null, ini, innovations = NULL,
                               Amat, rho_type, del = 1e-10, Bartlett = FALSE,
                               wDecay = FALSE, trunc_lev = TRUE,
                               sigvMethod = "factored",
                               winsorize_eps = FALSE) {
  p <- length(mdl_alt$phi)
  phi_cand <- theta[1:p]
  sigy_cand <- theta[p + 1]
  sigv_cand <- theta[p + 2]
  nu_cand <- theta[p + 3]
  stationary <- all(Mod(polyroot(c(1, -phi_cand))) > 1)
  if (!stationary || sigy_cand <= 0 || sigv_cand <= 0 || nu_cand <= 0) {
    return(9999999999999)
  }
  Tsize <- length(as.numeric(y))
  betasim_null <- c(phi_cand, sigy_cand, sigv_cand, nu_cand, rho_null)
  if (isTRUE(Bartlett)) {
    sN <- .simnull_lev_ged_Amat(betasim_null, rho_null, p, j, Tsize, N, ini,
                                 rho_type, del, TRUE, wDecay = wDecay,
                                 trunc_lev = trunc_lev,
                                 sigvMethod = sigvMethod,
                                 winsorize_eps = winsorize_eps,
                                 innovations = innovations)
  } else {
    sN <- .simnull_lev_ged(betasim_null, rho_null, p, j, Tsize, N, ini,
                            Amat, rho_type, del, wDecay = wDecay,
                            trunc_lev = trunc_lev,
                            sigvMethod = sigvMethod,
                            winsorize_eps = winsorize_eps,
                            innovations = innovations)
  }
  pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  return(pval)
}

# MMC p-value for Student-t (returns negative for minimization)
# S0 is kept fixed per Dufour (2006, eq 4.22)
.mmc_pval_t <- function(theta, y, j, N, s0, mdl_alt, nu_null, ini,
                        Amat, WAmat = FALSE, del = 1e-10, Bartlett = TRUE,
                        logNu = TRUE, wDecay = FALSE,
                        sigvMethod = "factored", winsorize_eps = 0,
                        direction = "two-sided",
                        innovations = NULL) {
  p <- length(mdl_alt$phi)
  Tsize <- length(as.numeric(y))
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (!stationary || theta[p + 1] <= 0 || theta[p + 2] <= 0) {
    return(9999999999999)
  }
  betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], nu_null)
  sN <- .simnull_t(betasim_null, nu_null, j, Tsize, N, ini,
                   Amat, del, WAmat, Bartlett, logNu,
                   wDecay = wDecay, sigvMethod = sigvMethod,
                   winsorize_eps = winsorize_eps,
                   direction = direction,
                   innovations = innovations)
  if (direction == "two-sided") {
    pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  } else {
    S_obs <- .signed_root(s0, mdl_alt$v, nu_null)
    pval <- -.pvalue_directional(S_obs, sN, direction)
  }
  return(pval)
}

# MMC p-value for GED (returns negative for minimization)
# S0 is kept fixed per Dufour (2006, eq 4.22)
.mmc_pval_ged <- function(theta, y, j, N, s0, mdl_alt, nu_null, ini,
                          Amat, WAmat = FALSE, del = 1e-10, Bartlett = TRUE,
                          wDecay = FALSE, sigvMethod = "factored",
                          winsorize_eps = 0,
                          innovations = NULL,
                          direction = "two-sided") {
  p <- length(mdl_alt$phi)
  Tsize <- length(as.numeric(y))
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (!stationary || theta[p + 1] <= 0 || theta[p + 2] <= 0) {
    return(9999999999999)
  }
  betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], nu_null)
  sN <- .simnull_ged(betasim_null, nu_null, j, Tsize, N, ini,
                     Amat, del, WAmat, Bartlett, wDecay = wDecay,
                     sigvMethod = sigvMethod,
                     winsorize_eps = winsorize_eps,
                     direction = direction,
                     innovations = innovations)
  if (direction == "two-sided") {
    pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  } else {
    S_obs <- .signed_root(s0, mdl_alt$v, nu_null)
    pval <- -.pvalue_directional(S_obs, sN, direction)
  }
  return(pval)
}

# --- Amat parsing and MMC optimizer dispatch ---

# Resolve Amat + Bartlett into a unified weighting specification.
# Amat takes precedence: "Weighted" -> HAC, <matrix> -> user-supplied, NULL -> check Bartlett.
# Bartlett = TRUE without Amat is backward-compat shorthand for Amat = "Weighted".
.resolve_weighting <- function(Amat, Bartlett) {
  if (!is.null(Amat)) {
    if (identical(Amat, "Weighted")) {
      return(list(Amat = "Weighted", use_hac = TRUE))
    } else if (is.matrix(Amat)) {
      return(list(Amat = Amat, use_hac = FALSE))
    } else {
      stop("Amat must be NULL, 'Weighted', or a numeric matrix.")
    }
  }
  if (isTRUE(Bartlett)) {
    return(list(Amat = "Weighted", use_hac = TRUE))
  }
  list(Amat = NULL, use_hac = FALSE)
}

# Parse Amat argument (shared by t, GED, and leverage tests)
# n_mom: number of moment conditions (p+3 for non-leverage heavy-tail, p+4 for leverage heavy-tail)
.parse_Amat <- function(Amat, p, n_mom = NULL) {
  if (is.null(n_mom)) n_mom <- p + 3
  WAmat <- FALSE
  if (is.null(Amat)) {
    Amat <- diag(n_mom)
  } else if (identical(Amat, "Weighted")) {
    WAmat <- TRUE
    Amat <- diag(n_mom)  # placeholder; overridden inside moment function
  } else if (!is.matrix(Amat) || !all(dim(Amat) == c(n_mom, n_mom))) {
    stop("Amat must be NULL, 'Weighted', or a (", n_mom, "x", n_mom, ") matrix.")
  }
  list(Amat = Amat, WAmat = WAmat)
}

# Generic MMC optimizer dispatch
.run_mmc_optimizer <- function(method, theta_0, fn, lower, upper,
                               threshold, maxit, ...) {
  if (method == "GenSA") {
    if (!requireNamespace("GenSA", quietly = TRUE)) {
      stop("Package 'GenSA' required. Install it or use method='pso'.")
    }
    if (is.null(maxit)) maxit <- 100
    out <- GenSA::GenSA(theta_0, fn, lower = lower, upper = upper,
                        control = list(threshold.stop = -threshold,
                                       max.call = maxit, verbose = TRUE),
                        ...)
  } else if (method == "pso") {
    if (!requireNamespace("pso", quietly = TRUE)) {
      stop("Package 'pso' required. Install it or use method='GenSA'.")
    }
    if (is.null(maxit)) maxit <- 100
    out <- pso::psoptim(theta_0, fn, lower = lower, upper = upper,
                        control = list(abstol = -threshold,
                                       maxf = maxit, trace = 1,
                                       REPORT = 1, trace.stats = TRUE),
                        ...)
  } else {
    stop("method must be 'pso' or 'GenSA'.")
  }
  return(out)
}

# --- GMM moment functions ---

# Pseudo-inverse via SVD
.pinv <- function(A, tol = .Machine$double.eps^0.5) {
  s <- svd(A)
  d <- s$d
  d[d > tol] <- 1 / d[d > tol]
  d[d <= tol] <- 0
  s$v %*% diag(d, nrow = length(d)) %*% t(s$u)
}

#' @keywords internal
LRT_moment_lev_svp <- function(y, mdl_out, Amat, rho_type, del = 1e-10) {
  LRT_moment_lev_svp_cpp(as.numeric(y), mdl_out, Amat, rho_type, del)
}

#' @keywords internal
LRT_moment_lev_svp_Amat <- function(y, mdl_out, rho_type, del = 1e-10,
                                     Bartlett = TRUE) {
  LRT_moment_lev_svp_Amat_cpp(as.numeric(y), mdl_out, rho_type, del,
                                Bartlett, pinv_fn = .pinv)
}

#' @keywords internal
LRT_moment_ar_Amat <- function(y, mdl_out, del = 1e-10, Bartlett = TRUE) {
  LRT_moment_ar_Amat_cpp(as.numeric(y), mdl_out, del, Bartlett, pinv_fn = .pinv)
}

#' @keywords internal
LRT_moment_t <- function(y, mdl_out, Amat, WAmat = FALSE, del = 1e-10,
                         Bartlett = TRUE) {
  LRT_moment_t_cpp(as.numeric(y), mdl_out, Amat, WAmat, del, Bartlett,
                    pinv_fn = .pinv)
}

#' @keywords internal
LRT_moment_ged <- function(y, mdl_out, Amat, WAmat = FALSE, del = 1e-10,
                           Bartlett = TRUE) {
  LRT_moment_ged_cpp(as.numeric(y), mdl_out, Amat, WAmat, del, Bartlett,
                      pinv_fn = .pinv)
}


# =========================================================================== #
# GMM moment functions for leverage + heavy-tail testing (p+4 conditions)
# =========================================================================== #

#' GMM moments for SVL(p)-Student-t with fixed A matrix (p+4 conditions)
#' @keywords internal
LRT_moment_lev_t <- function(y, mdl_out, Amat, rho_type, del = 1e-10) {
  LRT_moment_lev_t_cpp(as.numeric(y), mdl_out, Amat, rho_type, del)
}

#' GMM moments for SVL(p)-Student-t with HAC weighting (p+4 conditions)
#' @keywords internal
LRT_moment_lev_t_Amat <- function(y, mdl_out, rho_type, del = 1e-10,
                                   Bartlett = TRUE) {
  LRT_moment_lev_t_Amat_cpp(as.numeric(y), mdl_out, rho_type, del,
                              Bartlett, pinv_fn = .pinv)
}

#' GMM moments for SVL(p)-GED with fixed A matrix (p+4 conditions)
#' Uses exact GED leverage moment.
#' @keywords internal
LRT_moment_lev_ged <- function(y, mdl_out, Amat, rho_type, del = 1e-10) {
  gh <- .get_gh(200L)
  LRT_moment_lev_ged_cpp(as.numeric(y), mdl_out, Amat, rho_type, del,
                          gh_nodes = gh$nodes, gh_weights = gh$weights)
}

#' GMM moments for SVL(p)-GED with HAC weighting (p+4 conditions)
#' @keywords internal
LRT_moment_lev_ged_Amat <- function(y, mdl_out, rho_type, del = 1e-10,
                                     Bartlett = TRUE) {
  gh <- .get_gh(200L)
  LRT_moment_lev_ged_Amat_cpp(as.numeric(y), mdl_out, rho_type, del,
                                Bartlett, pinv_fn = .pinv,
                                gh_nodes = gh$nodes, gh_weights = gh$weights)
}



# =========================================================================== #
# Simulation-under-null helpers
# =========================================================================== #

# Simulate null distribution for AR test
.simnull_ar <- function(betasim_null, p_null, p_alt, j, Tsize, N, ini,
                        del = 1e-10, wDecay = FALSE, Bartlett = FALSE,
                        sigvMethod = "factored",
                        errorType = "Gaussian",
                        nu_fixed = NA_real_,
                        logNu = TRUE,
                        winsorize_eps = 0,
                        innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p_null)]
  sigy_sim <- betasim_null[p_null + 1]
  sigv_sim <- betasim_null[p_null + 2]
  n_total <- ini + Tsize

  # Pre-draw innovations if not supplied (errorType-specific eps distribution)
  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      eta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      eps = if (errorType == "Gaussian") {
        matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw)
      } else if (errorType == "Student-t") {
        matrix(rt(n_total * N_draw, df = nu_fixed), nrow = n_total, ncol = N_draw)
      } else {  # GED
        matrix(rged_cpp(n_total * N_draw, 0, 1, nu_fixed), nrow = n_total, ncol = N_draw)
      }
    )
  }

  # Helper closures for distribution-specific simulation, fitting, and stat computation
  sim_one <- function(eta_vec, eps_vec) {
    if (errorType == "Gaussian") {
      sim_from_innov_gaussian_cpp(phi = phi_sim, sigma_y = sigy_sim,
                                  sigma_v = sigv_sim,
                                  eta_vec = eta_vec, eps_vec = eps_vec,
                                  p = p_null, T_out = Tsize, burnin = ini)
    } else if (errorType == "Student-t") {
      sim_from_innov_t_cpp(phi = phi_sim, sigma_y = sigy_sim,
                           sigma_v = sigv_sim,
                           eta_vec = eta_vec, eps_vec = eps_vec,
                           p = p_null, T_out = Tsize, burnin = ini)
    } else {
      sim_from_innov_ged_cpp(phi = phi_sim, sigma_y = sigy_sim,
                             sigma_v = sigv_sim,
                             eta_vec = eta_vec, eps_vec = eps_vec,
                             p = p_null, T_out = Tsize, burnin = ini)
    }
  }
  fit_one <- function(u, p_use) {
    if (errorType == "Gaussian") {
      svp(u, p = p_use, J = j, leverage = FALSE, del = del, wDecay = wDecay,
          sigvMethod = sigvMethod)
    } else {
      svp(u, p = p_use, J = j, leverage = FALSE, errorType = errorType,
          del = del, wDecay = wDecay, logNu = logNu, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps)
    }
  }
  stat_bartlett <- function(u_vec, mdl_null_tmp, mdl_alt_tmp) {
    if (errorType == "Gaussian") {
      M_n <- LRT_moment_ar_Amat(u_vec, mdl_null_tmp, del = del, Bartlett = TRUE)
      M_a <- LRT_moment_ar_Amat(u_vec, mdl_alt_tmp, del = del, Bartlett = TRUE)
    } else {
      u_mat <- as.matrix(u_vec)
      wa_n <- .parse_Amat("Weighted", p_null)
      wa_a <- .parse_Amat("Weighted", p_alt)
      if (errorType == "Student-t") {
        M_n <- LRT_moment_t(u_mat, mdl_null_tmp, wa_n$Amat, wa_n$WAmat, del, TRUE)
        M_a <- LRT_moment_t(u_mat, mdl_alt_tmp, wa_a$Amat, wa_a$WAmat, del, TRUE)
      } else {  # GED
        M_n <- LRT_moment_ged(u_mat, mdl_null_tmp, wa_n$Amat, wa_n$WAmat, del, TRUE)
        M_a <- LRT_moment_ged(u_mat, mdl_alt_tmp, wa_a$Amat, wa_a$WAmat, del, TRUE)
      }
    }
    Tsize * (M_n - M_a)
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  while (xn <= N && b <= ncol(innovations$eta)) {
    u <- as.numeric(sim_one(innovations$eta[, b], innovations$eps[, b]))
    mdl_alt_tmp <- tryCatch(fit_one(u, p_alt), error = function(e) NULL)
    b <- b + 1
    if (is.null(mdl_alt_tmp) || isTRUE(mdl_alt_tmp$nonstationary_ind)) next
    if (isTRUE(Bartlett)) {
      mdl_null_tmp <- tryCatch(fit_one(u, p_null), error = function(e) NULL)
      if (is.null(mdl_null_tmp) || isTRUE(mdl_null_tmp$nonstationary_ind)) next
      stat <- tryCatch(stat_bartlett(u, mdl_null_tmp, mdl_alt_tmp),
                       error = function(e) NA_real_)
      if (is.na(stat)) next
      sN[xn] <- max(stat, 1e-10)
    } else {
      phi_extra <- mdl_alt_tmp$phi[(p_null + 1):p_alt]
      sN[xn] <- Tsize * sum(phi_extra^2)
    }
    xn <- xn + 1
  }
  # Fallback: draw fresh innovations if pre-drawn pool exhausted
  while (xn <= N) {
    eta_extra <- rnorm(n_total)
    eps_extra <- if (errorType == "Gaussian") rnorm(n_total)
                 else if (errorType == "Student-t") rt(n_total, df = nu_fixed)
                 else rged_cpp(n_total, 0, 1, nu_fixed)
    u <- as.numeric(sim_one(eta_extra, eps_extra))
    mdl_alt_tmp <- tryCatch(fit_one(u, p_alt), error = function(e) NULL)
    if (is.null(mdl_alt_tmp) || isTRUE(mdl_alt_tmp$nonstationary_ind)) next
    if (isTRUE(Bartlett)) {
      mdl_null_tmp <- tryCatch(fit_one(u, p_null), error = function(e) NULL)
      if (is.null(mdl_null_tmp) || isTRUE(mdl_null_tmp$nonstationary_ind)) next
      stat <- tryCatch(stat_bartlett(u, mdl_null_tmp, mdl_alt_tmp),
                       error = function(e) NA_real_)
      if (is.na(stat)) next
      sN[xn] <- max(stat, 1e-10)
    } else {
      phi_extra <- mdl_alt_tmp$phi[(p_null + 1):p_alt]
      sN[xn] <- Tsize * sum(phi_extra^2)
    }
    xn <- xn + 1
  }
  attr(sN, "innovations") <- innovations
  return(sN)
}

# Simulate null distribution for leverage test (identity Amat)
.simnull <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                     Amat, rho_type, del = 1e-10,
                     wDecay = FALSE, trunc_lev = TRUE,
                     sigvMethod = "factored",
                     innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  rho_sim <- betasim_null[p + 3]
  n_total <- ini + Tsize

  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      zeta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      aux  = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  max_iter <- ncol(innovations$zeta)
  while (xn <= N && b <= max_iter) {
    u <- as.numeric(sim_from_innov_gaussian_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim, rho = rho_sim,
      zeta_vec = innovations$zeta[, b], aux_vec = innovations$aux[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    b <- b + 1
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, rho_type = rho_type, del = del,
          trunc_lev = trunc_lev, wDecay = wDecay, sigvMethod = sigvMethod),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && abs(out_tmp$rho) <= 1 &&
        !isTRUE(out_tmp$nonstationary_ind)) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- Tsize * (LRT_moment_lev_svp(u_mat, out_null_tmp, Amat, rho_type, del) -
                           LRT_moment_lev_svp(u_mat, out_tmp, Amat, rho_type, del))
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  # Fallback
  while (xn <= N) {
    u <- as.numeric(sim_from_innov_gaussian_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim, rho = rho_sim,
      zeta_vec = rnorm(n_total), aux_vec = rnorm(n_total),
      p = p, T_out = Tsize, burnin = ini
    ))
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, rho_type = rho_type, del = del,
          trunc_lev = trunc_lev, wDecay = wDecay, sigvMethod = sigvMethod),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && abs(out_tmp$rho) <= 1 &&
        !isTRUE(out_tmp$nonstationary_ind)) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- Tsize * (LRT_moment_lev_svp(u_mat, out_null_tmp, Amat, rho_type, del) -
                           LRT_moment_lev_svp(u_mat, out_tmp, Amat, rho_type, del))
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  attr(sN, "innovations") <- innovations
  return(sN)
}

# Simulate null distribution for leverage test (Bartlett kernel Amat)
.simnull_Amat <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                          rho_type, del = 1e-10, Bartlett = TRUE,
                          wDecay = FALSE, trunc_lev = TRUE,
                          sigvMethod = "factored",
                          innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  rho_sim <- betasim_null[p + 3]
  n_total <- ini + Tsize

  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      zeta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      aux  = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  max_iter <- ncol(innovations$zeta)
  while (xn <= N && b <= max_iter) {
    u <- as.numeric(sim_from_innov_gaussian_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim, rho = rho_sim,
      zeta_vec = innovations$zeta[, b], aux_vec = innovations$aux[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    b <- b + 1
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, rho_type = rho_type, del = del,
          trunc_lev = trunc_lev, wDecay = wDecay, sigvMethod = sigvMethod),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && abs(out_tmp$rho) <= 1 &&
        !isTRUE(out_tmp$nonstationary_ind)) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- Tsize * (LRT_moment_lev_svp_Amat(u_mat, out_null_tmp, rho_type, del, Bartlett) -
                           LRT_moment_lev_svp_Amat(u_mat, out_tmp, rho_type, del, Bartlett))
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  while (xn <= N) {
    u <- as.numeric(sim_from_innov_gaussian_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim, rho = rho_sim,
      zeta_vec = rnorm(n_total), aux_vec = rnorm(n_total),
      p = p, T_out = Tsize, burnin = ini
    ))
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, rho_type = rho_type, del = del,
          trunc_lev = trunc_lev, wDecay = wDecay, sigvMethod = sigvMethod),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && abs(out_tmp$rho) <= 1 &&
        !isTRUE(out_tmp$nonstationary_ind)) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- Tsize * (LRT_moment_lev_svp_Amat(u_mat, out_null_tmp, rho_type, del, Bartlett) -
                           LRT_moment_lev_svp_Amat(u_mat, out_tmp, rho_type, del, Bartlett))
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  attr(sN, "innovations") <- innovations
  return(sN)
}

# Simulate null distribution for Student-t test
# Signed root: S = sign(theta_hat - theta_0) * sqrt(max(LR, 0))
.signed_root <- function(LR_T, theta_hat, theta_0) {
  sign(theta_hat - theta_0) * sqrt(pmax(LR_T, 0))
}

# Directional p-value
.pvalue_directional <- function(S_obs, S_sim, direction) {
  N <- length(S_sim)
  switch(direction,
    "less"    = (1 + sum(S_sim <= S_obs)) / (N + 1),
    "greater" = (1 + sum(S_sim >= S_obs)) / (N + 1)
  )
}

.simnull_t <- function(betasim_null, nu_null, j, Tsize, N, ini,
                       Amat, del = 1e-10, WAmat = FALSE,
                       Bartlett = TRUE, logNu = TRUE,
                       wDecay = FALSE, sigvMethod = "factored",
                       winsorize_eps = 0,
                       direction = "two-sided",
                       innovations = NULL) {
  p <- length(betasim_null) - 3
  phi_sim <- betasim_null[1:p]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  n_total <- ini + Tsize

  # Pre-draw innovations if not provided
  # When innovations is NULL: draw fresh (standard LMC or standalone MMC)
  # When innovations is provided: reuse (fixed-error MMC, Dufour 2006 eq 4.22)
  if (is.null(innovations)) {
    # Draw extra columns to handle estimation failures (need >= N valid)
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      eta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      eps = matrix(rt(n_total * N_draw, df = nu_sim), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  sign_vec <- numeric(N)
  xn <- 1
  b <- 1  # innovation index
  while (xn <= N && b <= ncol(innovations$eta)) {
    # Simulate from pre-drawn innovations
    u <- as.matrix(sim_from_innov_t_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      eta_vec = innovations$eta[, b], eps_vec = innovations$eps[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = j, errorType = "Student-t", del = del,
          logNu = logNu, wDecay = wDecay, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp)) {
      out_null_tmp <- out_tmp
      out_null_tmp$v <- nu_null
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_t(u, out_null_tmp, Amat, WAmat, del, Bartlett) -
                   LRT_moment_t(u, out_tmp, Amat, WAmat, del, Bartlett)),
        error = function(e) NULL
      )
      if (!is.null(sN_tmp) && !is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        sign_vec[xn] <- sign(out_tmp$v - nu_null)
        xn <- xn + 1
      }
    }
    b <- b + 1
  }
  # If we ran out of innovations, draw more (fallback)
  while (xn <= N) {
    eta_extra <- rnorm(n_total)
    eps_extra <- rt(n_total, df = nu_sim)
    u <- as.matrix(sim_from_innov_t_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      eta_vec = eta_extra, eps_vec = eps_extra,
      p = p, T_out = Tsize, burnin = ini
    ))
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = j, errorType = "Student-t", del = del,
          logNu = logNu, wDecay = wDecay, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp)) {
      out_null_tmp <- out_tmp
      out_null_tmp$v <- nu_null
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_t(u, out_null_tmp, Amat, WAmat, del, Bartlett) -
                   LRT_moment_t(u, out_tmp, Amat, WAmat, del, Bartlett)),
        error = function(e) NULL
      )
      if (!is.null(sN_tmp) && !is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        sign_vec[xn] <- sign(out_tmp$v - nu_null)
        xn <- xn + 1
      }
    }
  }

  out <- if (direction != "two-sided") sign_vec * sqrt(sN) else sN
  attr(out, "innovations") <- innovations
  return(out)
}

# Simulate null distribution for GED test
.simnull_ged <- function(betasim_null, nu_null, j, Tsize, N, ini,
                         Amat, del = 1e-10, WAmat = FALSE,
                         Bartlett = TRUE, wDecay = FALSE,
                         sigvMethod = "factored", winsorize_eps = 0,
                         direction = "two-sided",
                         innovations = NULL) {
  p <- length(betasim_null) - 3
  phi_sim <- betasim_null[1:p]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  n_total <- ini + Tsize

  # Pre-draw innovations if not provided
  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    a_ged <- exp(0.5 * (lgamma(1 / nu_sim) - lgamma(3 / nu_sim)))
    total_n <- n_total * N_draw
    x_gam_all <- rgamma(total_n, shape = 1 / nu_sim, rate = 1)
    eps_all <- sign(runif(total_n) - 0.5) * a_ged * x_gam_all^(1 / nu_sim)
    innovations <- list(
      eta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      eps = matrix(eps_all, nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  sign_vec <- numeric(N)
  xn <- 1
  b <- 1
  while (xn <= N && b <= ncol(innovations$eta)) {
    u <- as.matrix(sim_from_innov_ged_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      eta_vec = innovations$eta[, b], eps_vec = innovations$eps[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = j, errorType = "GED", del = del,
          wDecay = wDecay, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp)) {
      out_null_tmp <- out_tmp
      out_null_tmp$v <- nu_null
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_ged(u, out_null_tmp, Amat, WAmat, del, Bartlett) -
                   LRT_moment_ged(u, out_tmp, Amat, WAmat, del, Bartlett)),
        error = function(e) NULL
      )
      if (!is.null(sN_tmp) && !is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        sign_vec[xn] <- sign(out_tmp$v - nu_null)
        xn <- xn + 1
      }
    }
    b <- b + 1
  }
  # Fallback: draw more if needed
  while (xn <= N) {
    eta_extra <- rnorm(n_total)
    a_ged <- exp(0.5 * (lgamma(1 / nu_sim) - lgamma(3 / nu_sim)))
    x_gam <- rgamma(n_total, shape = 1 / nu_sim, rate = 1)
    eps_extra <- sign(runif(n_total) - 0.5) * a_ged * x_gam^(1 / nu_sim)
    u <- as.matrix(sim_from_innov_ged_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      eta_vec = eta_extra, eps_vec = eps_extra,
      p = p, T_out = Tsize, burnin = ini
    ))
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = j, errorType = "GED", del = del,
          wDecay = wDecay, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp)) {
      out_null_tmp <- out_tmp
      out_null_tmp$v <- nu_null
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_ged(u, out_null_tmp, Amat, WAmat, del, Bartlett) -
                   LRT_moment_ged(u, out_tmp, Amat, WAmat, del, Bartlett)),
        error = function(e) NULL
      )
      if (!is.null(sN_tmp) && !is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        sign_vec[xn] <- sign(out_tmp$v - nu_null)
        xn <- xn + 1
      }
    }
  }

  out <- if (direction != "two-sided") sign_vec * sqrt(sN) else sN
  attr(out, "innovations") <- innovations
  return(out)
}

# =========================================================================== #
# Null simulation helpers for leverage tests under heavy tails
# =========================================================================== #

.simnull_lev_t <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                            Amat, rho_type, del = 1e-10,
                            wDecay = FALSE, trunc_lev = TRUE,
                            logNu = FALSE, sigvMethod = "factored",
                            winsorize_eps = FALSE,
                            innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  rho_sim <- betasim_null[p + 4]
  n_total <- ini + Tsize

  # Student-t leverage: nu varies in nuisance → PIT for chi2
  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      zeta    = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      aux     = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      U_chi2  = matrix(runif(n_total * N_draw), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  max_iter <- ncol(innovations$zeta)
  while (xn <= N && b <= max_iter) {
    u <- as.numeric(sim_from_innov_t_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      nu = nu_sim, rho = rho_sim,
      zeta_vec = innovations$zeta[, b], aux_vec = innovations$aux[, b],
      U_chi2_vec = innovations$U_chi2[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    b <- b + 1
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, errorType = "Student-t",
          rho_type = rho_type, del = del, trunc_lev = trunc_lev,
          wDecay = wDecay, logNu = logNu, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && is.finite(out_tmp$v) &&
        !is.na(out_tmp$rho) && abs(out_tmp$rho) <= 1) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_lev_t(u_mat, out_null_tmp, Amat, rho_type, del) -
                   LRT_moment_lev_t(u_mat, out_tmp, Amat, rho_type, del)),
        error = function(e) NA
      )
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  if (xn <= N) warning("simnull_lev_t: only ", xn - 1, " of ", N, " valid stats.")
  attr(sN, "innovations") <- innovations
  return(sN)
}

.simnull_lev_t_Amat <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                                 rho_type, del = 1e-10, Bartlett = TRUE,
                                 wDecay = FALSE, trunc_lev = TRUE,
                                 logNu = FALSE, sigvMethod = "factored",
                                 winsorize_eps = FALSE,
                                 innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  rho_sim <- betasim_null[p + 4]
  n_total <- ini + Tsize

  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      zeta    = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      aux     = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      U_chi2  = matrix(runif(n_total * N_draw), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  max_iter <- ncol(innovations$zeta)
  while (xn <= N && b <= max_iter) {
    u <- as.numeric(sim_from_innov_t_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      nu = nu_sim, rho = rho_sim,
      zeta_vec = innovations$zeta[, b], aux_vec = innovations$aux[, b],
      U_chi2_vec = innovations$U_chi2[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    b <- b + 1
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, errorType = "Student-t",
          rho_type = rho_type, del = del, trunc_lev = trunc_lev,
          wDecay = wDecay, logNu = logNu, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && is.finite(out_tmp$v) &&
        !is.na(out_tmp$rho) && abs(out_tmp$rho) <= 1) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_lev_t_Amat(u_mat, out_null_tmp, rho_type, del, Bartlett) -
                   LRT_moment_lev_t_Amat(u_mat, out_tmp, rho_type, del, Bartlett)),
        error = function(e) NA
      )
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  if (xn <= N) warning("simnull_lev_t_Amat: only ", xn - 1, " of ", N, " valid stats.")
  attr(sN, "innovations") <- innovations
  return(sN)
}

.simnull_lev_ged <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                              Amat, rho_type, del = 1e-10,
                              wDecay = FALSE, trunc_lev = TRUE,
                              sigvMethod = "factored", winsorize_eps = FALSE,
                              innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  rho_sim <- betasim_null[p + 4]
  n_total <- ini + Tsize

  # GED leverage: Gaussian copula, just zeta+aux (no chi2 needed)
  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      zeta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      aux  = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  max_iter <- ncol(innovations$zeta)
  while (xn <= N && b <= max_iter) {
    u <- as.numeric(sim_from_innov_ged_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      nu = nu_sim, rho = rho_sim,
      zeta_vec = innovations$zeta[, b], aux_vec = innovations$aux[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    b <- b + 1
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, errorType = "GED",
          rho_type = rho_type, del = del, trunc_lev = trunc_lev,
          wDecay = wDecay, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && is.finite(out_tmp$v) &&
        !is.na(out_tmp$rho) && abs(out_tmp$rho) <= 1) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_lev_ged(u_mat, out_null_tmp, Amat, rho_type, del) -
                   LRT_moment_lev_ged(u_mat, out_tmp, Amat, rho_type, del)),
        error = function(e) NA
      )
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  if (xn <= N) warning("simnull_lev_ged: only ", xn - 1, " of ", N, " valid stats.")
  attr(sN, "innovations") <- innovations
  return(sN)
}

.simnull_lev_ged_Amat <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                                   rho_type, del = 1e-10, Bartlett = TRUE,
                                   wDecay = FALSE, trunc_lev = TRUE,
                                   sigvMethod = "factored", winsorize_eps = FALSE,
                                   innovations = NULL) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  rho_sim <- betasim_null[p + 4]
  n_total <- ini + Tsize

  if (is.null(innovations)) {
    N_draw <- ceiling(N * 1.5) + 10L
    innovations <- list(
      zeta = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw),
      aux  = matrix(rnorm(n_total * N_draw), nrow = n_total, ncol = N_draw)
    )
  }

  sN <- numeric(N)
  xn <- 1
  b <- 1
  max_iter <- ncol(innovations$zeta)
  while (xn <= N && b <= max_iter) {
    u <- as.numeric(sim_from_innov_ged_lev_cpp(
      phi = phi_sim, sigma_y = sigy_sim, sigma_v = sigv_sim,
      nu = nu_sim, rho = rho_sim,
      zeta_vec = innovations$zeta[, b], aux_vec = innovations$aux[, b],
      p = p, T_out = Tsize, burnin = ini
    ))
    b <- b + 1
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, errorType = "GED",
          rho_type = rho_type, del = del, trunc_lev = trunc_lev,
          wDecay = wDecay, sigvMethod = sigvMethod,
          winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && is.finite(out_tmp$v) &&
        !is.na(out_tmp$rho) && abs(out_tmp$rho) <= 1) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- tryCatch(
        Tsize * (LRT_moment_lev_ged_Amat(u_mat, out_null_tmp, rho_type, del, Bartlett) -
                   LRT_moment_lev_ged_Amat(u_mat, out_tmp, rho_type, del, Bartlett)),
        error = function(e) NA
      )
      if (!is.na(sN_tmp)) {
        sN[xn] <- max(sN_tmp, 1e-10)
        xn <- xn + 1
      }
    }
  }
  if (xn <= N) warning("simnull_lev_ged_Amat: only ", xn - 1, " of ", N, " valid stats.")
  attr(sN, "innovations") <- innovations
  return(sN)
}
