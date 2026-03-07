# =========================================================================== #
# Internal estimation helpers called by svp() in estim.R
# =========================================================================== #

# --- Gaussian SV(p) estimation (with optional leverage) ---
.svp_gaussian <- function(y, p, J, leverage, rho_type, del, trunc_lev, wDecay) {
  y <- as.numeric(y)
  if (length(y) < 2 * p + J) {
    stop("Time series too short for the given p and J.")
  }
  if (!rho_type %in% c("pearson", "kendall", "both", "none")) {
    stop("rho_type must be one of 'pearson', 'kendall', 'both', 'none'.")
  }
  if (!leverage) {
    rho_type <- "none"
  }
  if (rho_type == "none") {
    para <- svpCpp_nolev(y, p, J, del, wDecay)
    para$rho <- NA_real_
    para$gammatilde <- NA_real_
    para$theta <- c(para$phi, para$sigy, para$sigv)
  } else if (rho_type == "kendall") {
    para <- svpCpp(y, p, J, trunc_lev, del, rho_type = 2L, wDecay)
    para$rho_type <- rho_type
    para$theta <- c(para$phi, para$sigy, para$sigv, para$rho)
  } else if (rho_type == "pearson") {
    para <- svpCpp(y, p, J, trunc_lev, del, rho_type = 1L, wDecay)
    para$rho_type <- rho_type
    para$theta <- c(para$phi, para$sigy, para$sigv, para$rho)
  } else if (rho_type == "both") {
    para_ken <- svpCpp(y, p, J, trunc_lev, del, rho_type = 2L, wDecay)
    para <- svpCpp(y, p, J, trunc_lev, del, rho_type = 1L, wDecay)
    para$rho_kendall <- para_ken$rho
    para$rho_pearson <- para$rho
    para$rho_type <- rho_type
    para$theta <- c(para$phi, para$sigy, para$sigv, para$rho)
  }
  para$y <- y
  para$p <- p
  para$J <- J
  para$leverage <- leverage
  para$del <- del
  class(para) <- "svp"
  return(para)
}

# --- Student-t SV(p) estimation ---
.svp_t <- function(y, p, J, del, wDecay, logNu) {
  y <- as.matrix(as.numeric(y))
  N <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  # Estimate phi via W-ARMA-SV (distribution-free)
  para_base <- svpCpp_nolev(as.numeric(y), p, J, del, wDecay)
  phi_reg <- as.numeric(para_base$phi)
  # Compute theoretical lag-1 autocorrelation of AR(p) process
  rho_w1 <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = 1)[2])
  # Compute sample autocovariance at lag 1
  gam1 <- (1 / (N - 1)) * as.numeric(t(ys[2:N, , drop = FALSE]) %*% ys[1:(N - 1), , drop = FALSE])
  # Estimate sigma_eps^2 = gamma(0) - gamma(1)/rho_w(1)
  se2 <- stats::var(as.numeric(ly2)) - gam1 / rho_w1
  se2b <- as.numeric(se2) - psigamma(0.5, 1)
  se2b <- max(1e-6, se2b)
  # Estimate nu via root-finding on trigamma
  if (logNu) {
    f_log <- function(logx) psigamma(exp(logx) / 2, 1) - se2b
    lower <- log(2.01)
    upper <- log(500)
    if (f_log(lower) * f_log(upper) < 0) {
      nuh <- exp(stats::uniroot(f_log, interval = c(lower, upper), tol = 1e-6)$root)
    } else {
      warning("Root-finding for nu failed; returning NA.")
      nuh <- NA_real_
    }
  } else {
    f <- function(x) psigamma(x / 2, 1) - se2b
    lower <- 2.01
    upper <- 500
    if (f(lower) * f(upper) < 0) {
      nuh <- stats::uniroot(f, interval = c(lower, upper), tol = 1e-6, maxiter = 1000)$root
    } else {
      warning("Root-finding for nu failed; returning NA.")
      nuh <- NA_real_
    }
  }
  # Estimate sigv
  var_log_sq_t <- psigamma(0.5, 1) + psigamma(nuh / 2, 1)
  if (p == 1L) {
    # AD2021 formula for p=1 (more stable)
    sv_reg <- sqrt(abs((1 - phi_reg^2) * (stats::var(as.numeric(ly2)) - var_log_sq_t)))
  } else {
    # General formula: gamma(0) - sum(phi_j * gamma(j)) - sigma_eps^2
    gam_vec <- numeric(p)
    for (k in seq_len(p)) {
      gam_vec[k] <- (1 / (N - 1)) * as.numeric(
        t(ys[(k + 1):N, , drop = FALSE]) %*% ys[1:(N - k), , drop = FALSE])
    }
    sv_reg <- sqrt(abs(stats::var(as.numeric(ly2)) - sum(phi_reg * gam_vec) - var_log_sq_t))
  }
  # Estimate sigy
  mu_log_sq_t <- psigamma(0.5, 0) - psigamma(nuh / 2, 0) + log(nuh - 2)
  if (p == 1L) {
    sigy <- sqrt(exp(mean(log(y^2)) - mu_log_sq_t))
  } else {
    sigy <- sqrt(exp(mean(log(y^2 + del)) - mu_log_sq_t))
  }
  theta <- c(phi_reg, sigy, sv_reg, nuh)
  out <- list(mu = mu, phi = phi_reg, sigv = as.numeric(sv_reg),
              sigy = as.numeric(sigy), v = nuh, theta = theta,
              y = as.numeric(y), J = J, p = as.integer(p), del = del)
  class(out) <- "svp_t"
  return(out)
}

# --- GED SV(p) estimation ---
.svp_ged <- function(y, p, J, del, wDecay) {
  y <- as.matrix(as.numeric(y))
  N <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  # Estimate phi (distribution-free)
  para_base <- svpCpp_nolev(as.numeric(y), p, J, del, wDecay)
  phi_reg <- as.numeric(para_base$phi)
  # Compute theoretical lag-1 autocorrelation of AR(p) process
  rho_w1 <- as.numeric(stats::ARMAacf(ar = phi_reg, lag.max = 1)[2])
  # Compute sample autocovariance at lag 1
  gam1 <- (1 / (N - 1)) * as.numeric(t(ys[2:N, , drop = FALSE]) %*% ys[1:(N - 1), , drop = FALSE])
  # Estimate sigma_eps^2 = gamma(0) - gamma(1)/rho_w(1)
  se2 <- as.numeric(stats::var(as.numeric(ly2)) - gam1 / rho_w1)
  # Estimate nu
  f <- function(x) ((2 / x)^2) * psigamma(1 / x, 1) - se2
  m1 <- sum(abs(ly2)) / N
  m2 <- sqrt(sum((abs(ly2) - m1)^2) / N)
  x0 <- m1 / m2
  lower <- max(1e-6, x0 / 5)
  upper <- x0 * 5
  if (f(lower) * f(upper) < 0) {
    ged_nuh <- stats::uniroot(f, interval = c(lower, upper), tol = 1e-6, maxiter = 1000)$root
  } else {
    warning("Root-finding for nu failed; returning NA.")
    ged_nuh <- NA_real_
  }
  # Estimate sigv
  var_log_sq_ged <- ((2 / ged_nuh)^2) * psigamma(1 / ged_nuh, 1)
  if (p == 1L) {
    # AD2021 formula for p=1 (more stable)
    sv_reg <- sqrt(abs((1 - phi_reg^2) * (stats::var(as.numeric(ly2)) - var_log_sq_ged)))
  } else {
    # General formula: gamma(0) - sum(phi_j * gamma(j)) - sigma_eps^2
    gam_vec <- numeric(p)
    for (k in seq_len(p)) {
      gam_vec[k] <- (1 / (N - 1)) * as.numeric(
        t(ys[(k + 1):N, , drop = FALSE]) %*% ys[1:(N - k), , drop = FALSE])
    }
    sv_reg <- sqrt(abs(stats::var(as.numeric(ly2)) - sum(phi_reg * gam_vec) - var_log_sq_ged))
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
              y = as.numeric(y), J = J, p = as.integer(p), del = del)
  class(out) <- "svp_ged"
  return(out)
}





#' @noRd
rged <- function(n, mean = 0, sd = 1, nu = 2) {
  as.numeric(rged_cpp(as.integer(n), mean, sd, nu))
}




# --- Internal SE helpers ---

.svpSE_gaussian <- function(object, n_sim, alpha, burnin) {
  p <- object$p
  Tsize <- length(object$y)
  has_lev <- isTRUE(object$leverage) && !is.na(object$rho)
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
                   sigv = object$sigv, burnin = burnin)
    }
    out_tmp <- tryCatch(
      svp(u, p = p, J = object$J, leverage = has_lev,
          rho_type = rho_type, del = object$del, trunc_lev = TRUE),
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
  n_params <- length(object$phi) + 3
  betamat <- matrix(0, n_sim, n_params)
  betasim <- c(object$phi, object$sigy, object$sigv, object$v)
  xn <- 1
  while (xn <= n_sim) {
    u_out <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                     sigv = object$sigv, errorType = "Student-t",
                     nu = object$v, burnin = burnin)
    out_tmp <- tryCatch(
      svp(as.numeric(u_out), p = object$p, J = object$J,
          errorType = "Student-t", del = object$del, logNu = logNu),
      error = function(e) NULL
    )
    if (!is.null(out_tmp)) {
      betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv, out_tmp$v)
      xn <- xn + 1
    }
  }
  .compute_se_ci(betamat, betasim, n_sim, n_params, alpha)
}

.svpSE_ged <- function(object, n_sim, alpha, burnin) {
  Tsize <- length(object$y)
  n_params <- length(object$phi) + 3
  betamat <- matrix(0, n_sim, n_params)
  betasim <- c(object$phi, object$sigy, object$sigv, object$v)
  xn <- 1
  while (xn <= n_sim) {
    u_out <- sim_svp(Tsize, phi = object$phi, sigy = object$sigy,
                     sigv = object$sigv, errorType = "GED",
                     nu = object$v, burnin = burnin)
    out_tmp <- tryCatch(
      svp(as.numeric(u_out), p = object$p, J = object$J,
          errorType = "GED", del = object$del),
      error = function(e) NULL
    )
    if (!is.null(out_tmp)) {
      betamat[xn, ] <- c(out_tmp$phi, out_tmp$sigy, out_tmp$sigv, out_tmp$v)
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
                         del = 1e-10, wDecay = FALSE, Bartlett = FALSE) {
  Tsize <- length(y)
  phi_null <- theta[1:p_null]
  sigy_null <- theta[p_null + 1]
  sigv_null <- theta[p_null + 2]
  # Check stationarity
  stationary <- all(Mod(polyroot(c(1, -phi_null))) > 1)
  if (!stationary || sigy_null <= 0 || sigv_null <= 0) {
    return(9999999999999)
  }
  betasim_null <- c(phi_null, sigy_null, sigv_null)
  sN <- .simnull_ar(betasim_null, p_null, p_alt, j, Tsize, N, ini,
                    del, wDecay, Bartlett)
  pval <- -((N + 1 - sum(s0 >= sN)) / (N + 1))
  return(pval)
}

# MMC p-value for leverage test (identity Amat)
.mmc_pval_lev <- function(theta, y, j, N, mdl_alt, rho_null, ini,
                          Amat, rho_type, del = 1e-10) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  p <- length(mdl_alt$phi)
  out_lev_null <- list(phi = theta[1:p],
                       sigy = theta[p + 1],
                       sigv = theta[p + 2],
                       rho = rho_null,
                       gammatilde = mdl_alt$gammatilde)
  s0_tmp <- Tsize * (LRT_moment_lev_svp(y, out_lev_null, Amat, rho_type, del) -
                       LRT_moment_lev_svp(y, mdl_alt, Amat, rho_type, del))
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (!is.na(s0_tmp) && s0_tmp >= 0 && stationary) {
    betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], rho_null)
    sN <- .simnull(betasim_null, rho_null, p, j, Tsize, N, ini,
                   Amat, rho_type, del)
    pval <- -((N + 1 - sum(s0_tmp >= sN)) / (N + 1))
  } else {
    pval <- 9999999999999
  }
  return(pval)
}

# MMC p-value for leverage test (Bartlett kernel Amat)
.mmc_pval_lev_Amat <- function(theta, y, j, N, mdl_alt, rho_null, ini,
                               rho_type, del = 1e-10, Bartlett = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  p <- length(mdl_alt$phi)
  out_lev_null <- list(phi = theta[1:p],
                       sigy = theta[p + 1],
                       sigv = theta[p + 2],
                       rho = rho_null,
                       gammatilde = mdl_alt$gammatilde)
  s0_tmp <- Tsize * (LRT_moment_lev_svp_Amat(y, out_lev_null, rho_type, del, Bartlett) -
                       LRT_moment_lev_svp_Amat(y, mdl_alt, rho_type, del, Bartlett))
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (!is.na(s0_tmp) && s0_tmp >= 0 && stationary) {
    betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], rho_null)
    sN <- .simnull_Amat(betasim_null, rho_null, p, j, Tsize, N, ini,
                        rho_type, del, Bartlett)
    pval <- -((N + 1 - sum(s0_tmp >= sN)) / (N + 1))
  } else {
    pval <- 9999999999999
  }
  return(pval)
}

# MMC p-value for Student-t (returns negative for minimization)
.mmc_pval_t <- function(theta, y, j, N, mdl_alt, nu_null, ini,
                        Amat, del = 1e-10, Bartlett = TRUE,
                        logNu = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  p <- length(mdl_alt$phi)
  out_lev_null <- list(phi = theta[1:p], sigy = theta[p + 1],
                       sigv = theta[p + 2], v = nu_null)
  s0_tmp <- tryCatch(
    Tsize * (LRT_moment_t(y, out_lev_null, Amat, FALSE, del, Bartlett) -
               LRT_moment_t(y, mdl_alt, Amat, FALSE, del, Bartlett)),
    error = function(e) NA
  )
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (is.finite(s0_tmp) && s0_tmp >= 0 && stationary) {
    betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], nu_null)
    sN <- .simnull_t(betasim_null, nu_null, j, Tsize, N, ini,
                     Amat, del, FALSE, Bartlett, logNu)
    pval <- -((N + 1 - sum(s0_tmp >= sN)) / (N + 1))
  } else {
    pval <- 9999999999999
  }
  return(pval)
}

# MMC p-value for GED (returns negative for minimization)
.mmc_pval_ged <- function(theta, y, j, N, mdl_alt, nu_null, ini,
                          Amat, del = 1e-10, Bartlett = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  p <- length(mdl_alt$phi)
  out_lev_null <- list(phi = theta[1:p], sigy = theta[p + 1],
                       sigv = theta[p + 2], v = nu_null)
  s0_tmp <- tryCatch(
    Tsize * (LRT_moment_ged(y, out_lev_null, Amat, FALSE, del, Bartlett) -
               LRT_moment_ged(y, mdl_alt, Amat, FALSE, del, Bartlett)),
    error = function(e) NA
  )
  stationary <- all(Mod(polyroot(c(1, -theta[1:p]))) > 1)
  if (is.finite(s0_tmp) && s0_tmp >= 0 && stationary) {
    betasim_null <- c(theta[1:p], theta[p + 1], theta[p + 2], nu_null)
    sN <- .simnull_ged(betasim_null, nu_null, j, Tsize, N, ini,
                       Amat, del, FALSE, Bartlett)
    pval <- -((N + 1 - sum(s0_tmp >= sN)) / (N + 1))
  } else {
    pval <- 9999999999999
  }
  return(pval)
}

# --- Amat parsing and MMC optimizer dispatch ---

# Parse Amat argument (shared by t and GED tests)
.parse_Amat <- function(Amat, p) {
  WAmat <- FALSE
  if (is.null(Amat)) {
    Amat <- diag(p + 3)
  } else if (identical(Amat, "Weighted")) {
    WAmat <- TRUE
    Amat <- diag(p + 3)  # placeholder; overridden inside moment function
  } else if (!is.matrix(Amat) || !all(dim(Amat) == c(p + 3, p + 3))) {
    stop("Amat must be NULL, 'Weighted', or a (p+3)x(p+3) matrix.")
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
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  phi <- as.numeric(mdl_out$phi)
  p <- length(phi)
  gam0 <- acov_g(ys, 0)
  gamk <- numeric(2 * p)
  for (xp in seq_len(2 * p)) {
    gamk[xp] <- acov_g(ys, xp)
  }
  gamtmp <- 0
  if (p >= 2) {
    for (xp in 2:p) {
      gamtmp <- gamtmp + phi[xp] * (gamk[xp - 1] + gamk[xp])
    }
  }
  yabs <- abs(y)
  yn <- y
  muu <- mean(yn[1:(Tsize - 1), ])
  mua <- mean(yabs[2:Tsize, ])
  if (rho_type == "kendall") {
    EH_kentmp <- kendall_corr(yabs[2:Tsize, ] - mua, yn[1:(Tsize - 1), ] - muu)
    EH <- EH_kentmp * sqrt(stats::var(yabs[2:Tsize, ]) * stats::var(yn[1:(Tsize - 1), ]))
  } else {
    EH <- (t(yabs[2:Tsize, , drop = FALSE] - mua) %*%
             (yn[1:(Tsize - 1), , drop = FALSE] - muu)) / (Tsize - 2)
  }
  m1 <- mu + 1.2704 - log(mdl_out$sigy^2)
  m2 <- gam0 + gamk[1] - ((pi^2) / 2) -
    (1 / (1 - phi[1])) * (gamtmp + mdl_out$sigv^2)
  mk <- numeric(0)
  for (xp in (p + 1):(2 * p)) {
    mk <- c(mk, gamk[xp] - sum(phi * gamk[xp - (1:p)]))
  }
  m4 <- mdl_out$rho - ((as.numeric(EH) * sqrt(2 * pi)) /
    (mdl_out$sigv * (mdl_out$sigy^2))) * exp(-0.25 * mdl_out$gammatilde)
  g <- as.matrix(c(m1, m2, mk, m4))
  M_lev <- as.numeric(t(g) %*% Amat %*% g)
  return(M_lev)
}

#' @keywords internal
LRT_moment_lev_svp_Amat <- function(y, mdl_out, rho_type, del = 1e-10,
                                     Bartlett = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  phi <- as.numeric(mdl_out$phi)
  p <- length(phi)
  g_t <- matrix(0, Tsize, p + 3)
  gam0 <- acov_g(ys, 0)
  gam0_t <- ys * ys
  gamk <- numeric(2 * p)
  gamk_t <- matrix(0, Tsize, 2 * p)
  for (xp in seq_len(2 * p)) {
    gamk[xp] <- acov_g(ys, xp)
    gamk_t[(xp + 1):Tsize, xp] <- ys[(xp + 1):Tsize, ] * ys[1:(Tsize - xp), ]
  }
  gamtmp <- 0
  gamtmp_t <- matrix(0, Tsize, 1)
  if (p >= 2) {
    for (xp in 2:p) {
      gamtmp <- gamtmp + phi[xp] * (gamk[xp - 1] + gamk[xp])
      gamtmp_t <- gamtmp_t + phi[xp] * (gamk_t[, xp - 1] + gamk_t[, xp])
    }
  }
  yabs <- abs(y)
  yn <- y
  muu <- mean(yn[1:(Tsize - 1), ])
  mua <- mean(yabs[2:Tsize, ])
  if (rho_type == "kendall") {
    EH_kentmp <- kendall_corr(yabs[2:Tsize, ] - mua, yn[1:(Tsize - 1), ] - muu)
    EH <- EH_kentmp * sqrt(stats::var(yabs[2:Tsize, ]) * stats::var(yn[1:(Tsize - 1), ]))
    g_t[2:Tsize, ncol(g_t)] <- EH_kentmp *
      sqrt((yabs[2:Tsize, ]^2) * (yn[1:(Tsize - 1), ]^2))
  } else {
    EH <- (t(yabs[2:Tsize, , drop = FALSE] - mua) %*%
             (yn[1:(Tsize - 1), , drop = FALSE] - muu)) / (Tsize - 2)
    g_t[2:Tsize, ncol(g_t)] <- (yabs[2:Tsize, , drop = FALSE]) *
      (yn[1:(Tsize - 1), , drop = FALSE])
  }
  m1 <- mu + 1.2704 - log(mdl_out$sigy^2)
  g_t[, 1] <- ly2 + 1.2704 - log(mdl_out$sigy^2)
  m2 <- gam0 + gamk[1] - ((pi^2) / 2) -
    (1 / (1 - phi[1])) * (gamtmp + mdl_out$sigv^2)
  g_t[, 2] <- gam0_t + gamk_t[, 1] - ((pi^2) / 2) -
    (1 / (1 - phi[1])) * (gamtmp_t + mdl_out$sigv^2)
  mk <- numeric(0)
  mk_t <- matrix(0, Tsize, 0)
  for (xp in (p + 1):(2 * p)) {
    mk <- c(mk, gamk[xp] - sum(phi * gamk[xp - (1:p)]))
    mk_t <- cbind(mk_t, gamk_t[, xp] -
      rowSums(gamk_t[, xp - (1:p), drop = FALSE] %*% as.matrix(phi)))
  }
  g_t[, 3:(ncol(g_t) - 1)] <- mk_t
  m4 <- mdl_out$rho - ((as.numeric(EH) * sqrt(2 * pi)) /
    (mdl_out$sigv * (mdl_out$sigy^2))) * exp(-0.25 * mdl_out$gammatilde)
  g_t[2:Tsize, ncol(g_t)] <- mdl_out$rho -
    ((g_t[2:Tsize, ncol(g_t)] * sqrt(2 * pi)) /
       (mdl_out$sigv * (mdl_out$sigy^2))) * exp(-0.25 * mdl_out$gammatilde)
  g_t <- g_t[(p + 1):nrow(g_t), , drop = FALSE]
  g <- as.matrix(c(m1, m2, mk, m4))
  Gam0 <- (1 / (Tsize - p)) * (t(g_t) %*% g_t)
  if (isTRUE(Bartlett)) {
    KT <- floor(Tsize^(1 / 3))
    W <- Gam0
    for (k in seq_len(KT)) {
      weight <- 1 - (k / (KT + 1))
      Gamk <- t(g_t[(k + 1):nrow(g_t), , drop = FALSE]) %*%
        g_t[1:(nrow(g_t) - k), , drop = FALSE] / (Tsize - 3)
      W <- W + weight * (Gamk + t(Gamk))
    }
    Amat <- .pinv(W)
  } else {
    Amat <- .pinv(Gam0)
  }
  M <- as.numeric(t(g) %*% Amat %*% g)
  return(M)
}

#' @keywords internal
LRT_moment_ar_Amat <- function(y, mdl_out, del = 1e-10, Bartlett = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  phi <- as.numeric(mdl_out$phi)
  p <- length(phi)
  # Number of moment conditions: 1 (mean) + 1 (variance) + p (autocovariance)
  n_mom <- p + 2
  g_t <- matrix(0, Tsize, n_mom)
  gam0 <- acov_g(ys, 0)
  gam0_t <- ys * ys
  gamk <- numeric(2 * p)
  gamk_t <- matrix(0, Tsize, 2 * p)
  for (xp in seq_len(2 * p)) {
    gamk[xp] <- acov_g(ys, xp)
    gamk_t[(xp + 1):Tsize, xp] <- ys[(xp + 1):Tsize, ] * ys[1:(Tsize - xp), ]
  }
  gamtmp <- 0
  gamtmp_t <- matrix(0, Tsize, 1)
  if (p >= 2) {
    for (xp in 2:p) {
      gamtmp <- gamtmp + phi[xp] * (gamk[xp - 1] + gamk[xp])
      gamtmp_t <- gamtmp_t + phi[xp] * (gamk_t[, xp - 1] + gamk_t[, xp])
    }
  }
  # Moment 1: mean condition
  m1 <- mu + 1.2704 - log(mdl_out$sigy^2)
  g_t[, 1] <- ly2 + 1.2704 - log(mdl_out$sigy^2)
  # Moment 2: variance condition
  m2 <- gam0 + gamk[1] - ((pi^2) / 2) -
    (1 / (1 - phi[1])) * (gamtmp + mdl_out$sigv^2)
  g_t[, 2] <- gam0_t + gamk_t[, 1] - ((pi^2) / 2) -
    (1 / (1 - phi[1])) * (gamtmp_t + mdl_out$sigv^2)
  # Moments 3+: autocovariance conditions
  col_idx <- 3
  for (xp in (p + 1):(2 * p)) {
    g_t[, col_idx] <- gamk_t[, xp] -
      rowSums(gamk_t[, xp - (1:p), drop = FALSE] %*% as.matrix(phi))
    col_idx <- col_idx + 1
  }
  mk <- numeric(0)
  for (xp in (p + 1):(2 * p)) {
    mk <- c(mk, gamk[xp] - sum(phi * gamk[xp - (1:p)]))
  }
  g <- as.matrix(c(m1, m2, mk))
  # Trim initial rows
  g_t <- g_t[(p + 1):nrow(g_t), , drop = FALSE]
  T_eff <- nrow(g_t)
  # HAC estimation
  Gam0 <- (1 / T_eff) * (t(g_t) %*% g_t)
  if (isTRUE(Bartlett)) {
    KT <- floor(Tsize^(1 / 3))
    W <- Gam0
    for (k in seq_len(KT)) {
      weight <- 1 - (k / (KT + 1))
      Gamk <- t(g_t[(k + 1):T_eff, , drop = FALSE]) %*%
        g_t[1:(T_eff - k), , drop = FALSE] / (T_eff - 1)
      W <- W + weight * (Gamk + t(Gamk))
    }
    Amat <- .pinv(W)
  } else {
    Amat <- .pinv(Gam0)
  }
  M <- as.numeric(t(g) %*% Amat %*% g)
  return(M)
}

#' @keywords internal
LRT_moment_t <- function(y, mdl_out, Amat, WAmat = FALSE, del = 1e-10,
                         Bartlett = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  phi <- as.numeric(mdl_out$phi)
  p <- length(phi)
  n_mom <- p + 3
  g_t <- matrix(0, Tsize, n_mom)
  gam0 <- acov_g(ys, 0)
  gam0_t <- ys * ys
  gamk <- numeric(2 * p)
  gamk_t <- matrix(0, Tsize, 2 * p)
  for (xp in seq_len(2 * p)) {
    gamk[xp] <- acov_g(ys, xp)
    gamk_t[(xp + 1):Tsize, xp] <- ys[(xp + 1):Tsize, ] * ys[1:(Tsize - xp), ]
  }
  # Compute theoretical lag-1 autocorrelation
  rho_w1 <- as.numeric(stats::ARMAacf(ar = phi, lag.max = 1)[2])
  # Distribution-specific constants
  mu_log_sq_t <- psigamma(0.5, 0) - psigamma(mdl_out$v / 2, 0) + log(mdl_out$v - 2)
  var_log_sq_t <- psigamma(0.5, 1) + psigamma(mdl_out$v / 2, 1)
  # g1: mean condition
  m1 <- mu - mu_log_sq_t - log(mdl_out$sigy^2)
  # g2: sigma_v^2 condition: gamma(0) - sum(phi_j * gamma(j)) - sigma_eps^2 - sigma_v^2
  m2 <- gam0 - sum(phi * gamk[1:p]) - var_log_sq_t - (mdl_out$sigv^2)
  # g3..g_{p+2}: YW overidentifying restrictions at lags p+1,...,2p
  mk <- numeric(p)
  for (xp in seq_len(p)) {
    lag_idx <- p + xp
    mk[xp] <- gamk[lag_idx] - sum(phi * gamk[lag_idx - (1:p)])
  }
  # g_{p+3}: profiling condition: sigma_eps^2(nu) - gamma(0) + gamma(1)/rho_w(1)
  m_prof <- var_log_sq_t - gam0 + gamk[1] / rho_w1
  g <- as.matrix(c(m1, m2, mk, m_prof))
  if (isTRUE(WAmat)) {
    g_t[, 1] <- ly2 - mu_log_sq_t - log(mdl_out$sigy^2)
    g_t[, 2] <- gam0_t - rowSums(gamk_t[, 1:p, drop = FALSE] %*% diag(phi, nrow = p)) -
      var_log_sq_t - (mdl_out$sigv^2)
    for (xp in seq_len(p)) {
      lag_idx <- p + xp
      g_t[, 2 + xp] <- gamk_t[, lag_idx] -
        rowSums(gamk_t[, lag_idx - (1:p), drop = FALSE] %*% diag(phi, nrow = p))
    }
    g_t[, n_mom] <- var_log_sq_t - gam0_t + gamk_t[, 1] / rho_w1
    g_t <- g_t[(p + 1):nrow(g_t), , drop = FALSE]
    T_eff <- nrow(g_t)
    Gam0 <- (1 / T_eff) * (t(g_t) %*% g_t)
    if (isTRUE(Bartlett)) {
      KT <- floor(Tsize^(1 / 3))
      W <- Gam0
      for (k in seq_len(KT)) {
        weight <- 1 - (k / (KT + 1))
        Gamk <- t(g_t[(k + 1):T_eff, , drop = FALSE]) %*%
          g_t[1:(T_eff - k), , drop = FALSE] / (T_eff - 1)
        W <- W + weight * (Gamk + t(Gamk))
      }
      Amat <- .pinv(W)
    } else {
      Amat <- .pinv(Gam0)
    }
  }
  M <- as.numeric(t(g) %*% Amat %*% g)
  return(M)
}

#' @keywords internal
LRT_moment_ged <- function(y, mdl_out, Amat, WAmat = FALSE, del = 1e-10,
                           Bartlett = TRUE) {
  y <- as.matrix(as.numeric(y))
  Tsize <- nrow(y)
  ly2 <- log(y^2 + del)
  mu <- mean(ly2)
  ys <- ly2 - mu
  phi <- as.numeric(mdl_out$phi)
  p <- length(phi)
  n_mom <- p + 3
  g_t <- matrix(0, Tsize, n_mom)
  gam0 <- acov_g(ys, 0)
  gam0_t <- ys * ys
  gamk <- numeric(2 * p)
  gamk_t <- matrix(0, Tsize, 2 * p)
  for (xp in seq_len(2 * p)) {
    gamk[xp] <- acov_g(ys, xp)
    gamk_t[(xp + 1):Tsize, xp] <- ys[(xp + 1):Tsize, ] * ys[1:(Tsize - xp), ]
  }
  # Compute theoretical lag-1 autocorrelation
  rho_w1 <- as.numeric(stats::ARMAacf(ar = phi, lag.max = 1)[2])
  # Distribution-specific constants
  mu_log_sq_ged <- (2 / mdl_out$v) * psigamma(1 / mdl_out$v, 0) +
    log(gamma(1 / mdl_out$v)) - log(gamma(3 / mdl_out$v))
  var_log_sq_ged <- ((2 / mdl_out$v)^2) * psigamma(1 / mdl_out$v, 1)
  # g1: mean condition
  m1 <- mu - mu_log_sq_ged - log(mdl_out$sigy^2)
  # g2: sigma_v^2 condition
  m2 <- gam0 - sum(phi * gamk[1:p]) - var_log_sq_ged - (mdl_out$sigv^2)
  # g3..g_{p+2}: YW overidentifying restrictions
  mk <- numeric(p)
  for (xp in seq_len(p)) {
    lag_idx <- p + xp
    mk[xp] <- gamk[lag_idx] - sum(phi * gamk[lag_idx - (1:p)])
  }
  # g_{p+3}: profiling condition
  m_prof <- var_log_sq_ged - gam0 + gamk[1] / rho_w1
  g <- as.matrix(c(m1, m2, mk, m_prof))
  if (isTRUE(WAmat)) {
    g_t[, 1] <- ly2 - mu_log_sq_ged - log(mdl_out$sigy^2)
    g_t[, 2] <- gam0_t - rowSums(gamk_t[, 1:p, drop = FALSE] %*% diag(phi, nrow = p)) -
      var_log_sq_ged - (mdl_out$sigv^2)
    for (xp in seq_len(p)) {
      lag_idx <- p + xp
      g_t[, 2 + xp] <- gamk_t[, lag_idx] -
        rowSums(gamk_t[, lag_idx - (1:p), drop = FALSE] %*% diag(phi, nrow = p))
    }
    g_t[, n_mom] <- var_log_sq_ged - gam0_t + gamk_t[, 1] / rho_w1
    g_t <- g_t[(p + 1):nrow(g_t), , drop = FALSE]
    T_eff <- nrow(g_t)
    Gam0 <- (1 / T_eff) * (t(g_t) %*% g_t)
    if (isTRUE(Bartlett)) {
      KT <- floor(Tsize^(1 / 3))
      W <- Gam0
      for (k in seq_len(KT)) {
        weight <- 1 - (k / (KT + 1))
        Gamk <- t(g_t[(k + 1):T_eff, , drop = FALSE]) %*%
          g_t[1:(T_eff - k), , drop = FALSE] / (T_eff - 1)
        W <- W + weight * (Gamk + t(Gamk))
      }
      Amat <- .pinv(W)
    } else {
      Amat <- .pinv(Gam0)
    }
  }
  M <- as.numeric(t(g) %*% Amat %*% g)
  return(M)
}





# =========================================================================== #
# Simulation-under-null helpers
# =========================================================================== #

# Simulate null distribution for AR test
.simnull_ar <- function(betasim_null, p_null, p_alt, j, Tsize, N, ini,
                        del = 1e-10, wDecay = FALSE, Bartlett = FALSE) {
  phi_sim <- betasim_null[seq_len(p_null)]
  sigy_sim <- betasim_null[p_null + 1]
  sigv_sim <- betasim_null[p_null + 2]
  sN <- numeric(N)
  xn <- 1
  while (xn <= N) {
    u <- sim_svp(Tsize, phi = phi_sim, sigy = sigy_sim,
                 sigv = sigv_sim, burnin = ini)
    mdl_alt_tmp <- tryCatch(
      svp(u, p = p_alt, J = j, leverage = FALSE, del = del, wDecay = wDecay),
      error = function(e) NULL
    )
    if (is.null(mdl_alt_tmp) || isTRUE(mdl_alt_tmp$nonstationary_ind)) next
    if (isTRUE(Bartlett)) {
      mdl_null_tmp <- tryCatch(
        svp(u, p = p_null, J = j, leverage = FALSE, del = del, wDecay = wDecay),
        error = function(e) NULL
      )
      if (is.null(mdl_null_tmp) || isTRUE(mdl_null_tmp$nonstationary_ind)) next
      M_null <- LRT_moment_ar_Amat(u, mdl_null_tmp, del = del, Bartlett = TRUE)
      M_alt  <- LRT_moment_ar_Amat(u, mdl_alt_tmp, del = del, Bartlett = TRUE)
      stat <- Tsize * (M_null - M_alt)
      if (stat < 0) next  # discard negative test stats
      sN[xn] <- stat
    } else {
      phi_extra <- mdl_alt_tmp$phi[(p_null + 1):p_alt]
      sN[xn] <- Tsize * sum(phi_extra^2)
    }
    xn <- xn + 1
  }
  return(sN)
}

# Simulate null distribution for leverage test (identity Amat)
.simnull <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                     Amat, rho_type, del = 1e-10) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  rho_sim <- betasim_null[p + 3]
  sN <- numeric(N)
  xn <- 1
  while (xn <= N) {
    u_out <- sim_svp(Tsize, phi = phi_sim, sigy = sigy_sim,
                     sigv = sigv_sim, leverage = TRUE,
                     rho = rho_sim, burnin = ini)
    u <- u_out$y
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, rho_type = rho_type, del = del),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && abs(out_tmp$rho) <= 1 &&
        !isTRUE(out_tmp$nonstationary_ind)) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- Tsize * (LRT_moment_lev_svp(u_mat, out_null_tmp, Amat, rho_type, del) -
                           LRT_moment_lev_svp(u_mat, out_tmp, Amat, rho_type, del))
      if (!is.na(sN_tmp) && sN_tmp >= 0) {
        sN[xn] <- sN_tmp
        xn <- xn + 1
      }
    }
  }
  return(sN)
}

# Simulate null distribution for leverage test (Bartlett kernel Amat)
.simnull_Amat <- function(betasim_null, rho_null, p, j, Tsize, N, ini,
                          rho_type, del = 1e-10, Bartlett = TRUE) {
  phi_sim <- betasim_null[seq_len(p)]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  rho_sim <- betasim_null[p + 3]
  sN <- numeric(N)
  xn <- 1
  while (xn <= N) {
    u_out <- sim_svp(Tsize, phi = phi_sim, sigy = sigy_sim,
                     sigv = sigv_sim, leverage = TRUE,
                     rho = rho_sim, burnin = ini)
    u <- u_out$y
    out_tmp <- tryCatch(
      svp(u, p, j, leverage = TRUE, rho_type = rho_type, del = del),
      error = function(e) NULL
    )
    if (!is.null(out_tmp) && abs(out_tmp$rho) <= 1 &&
        !isTRUE(out_tmp$nonstationary_ind)) {
      out_null_tmp <- out_tmp
      out_null_tmp$rho <- rho_null
      u_mat <- as.matrix(u)
      sN_tmp <- Tsize * (LRT_moment_lev_svp_Amat(u_mat, out_null_tmp, rho_type, del, Bartlett) -
                           LRT_moment_lev_svp_Amat(u_mat, out_tmp, rho_type, del, Bartlett))
      if (!is.na(sN_tmp) && sN_tmp >= 0) {
        sN[xn] <- sN_tmp
        xn <- xn + 1
      }
    }
  }
  return(sN)
}

# Simulate null distribution for Student-t test
.simnull_t <- function(betasim_null, nu_null, j, Tsize, N, ini,
                       Amat, del = 1e-10, WAmat = FALSE,
                       Bartlett = TRUE, logNu = TRUE) {
  p <- length(betasim_null) - 3
  phi_sim <- betasim_null[1:p]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  sN <- numeric(N)
  xn <- 1
  while (xn <= N) {
    u <- as.matrix(sim_svp(Tsize, phi = phi_sim, sigy = sigy_sim,
                           sigv = sigv_sim, errorType = "Student-t",
                           nu = nu_sim, burnin = ini))
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = j, errorType = "Student-t", del = del, logNu = logNu),
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
      if (!is.null(sN_tmp) && !is.na(sN_tmp) && sN_tmp >= 0) {
        sN[xn] <- sN_tmp
        xn <- xn + 1
      }
    }
  }
  return(sN)
}

# Simulate null distribution for GED test
.simnull_ged <- function(betasim_null, nu_null, j, Tsize, N, ini,
                         Amat, del = 1e-10, WAmat = FALSE,
                         Bartlett = TRUE) {
  p <- length(betasim_null) - 3
  phi_sim <- betasim_null[1:p]
  sigy_sim <- betasim_null[p + 1]
  sigv_sim <- betasim_null[p + 2]
  nu_sim <- betasim_null[p + 3]
  sN <- numeric(N)
  xn <- 1
  while (xn <= N) {
    u <- as.matrix(sim_svp(Tsize, phi = phi_sim, sigy = sigy_sim,
                           sigv = sigv_sim, errorType = "GED",
                           nu = nu_sim, burnin = ini))
    out_tmp <- tryCatch(
      svp(as.numeric(u), p = p, J = j, errorType = "GED", del = del),
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
      if (!is.null(sN_tmp) && !is.na(sN_tmp) && sN_tmp >= 0) {
        sN[xn] <- sN_tmp
        xn <- xn + 1
      }
    }
  }
  return(sN)
}
