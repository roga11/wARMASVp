#' Multi-Step Ahead Volatility Forecast
#'
#' Applies Kalman filtering/smoothing to an estimated SV(p) model and produces
#' multi-step ahead volatility forecasts with uncertainty quantification.
#'
#' @param object An \code{"svp"}, \code{"svp_t"}, or \code{"svp_ged"} object
#'   from \code{\link{svp}}.
#' @param H Integer. Maximum forecast horizon. Default 1.
#' @param output Character. Primary output scale: \code{"log-variance"} (default,
#'   native log-volatility w_h), \code{"variance"} (conditional variance
#'   \eqn{\sigma^2_{T+h|T}}), or \code{"volatility"} (conditional std dev
#'   \eqn{\sigma_{T+h|T}}). All three are always computed and stored; this
#'   controls which is used by \code{print} and \code{plot} methods.
#' @param filter_method Character. Filter method: \code{"corrected"} (default),
#'   \code{"mixture"} (GMKF), or \code{"particle"} (BPF).
#' @param mixture Character. Source of the GMKF measurement mixture when
#'   \code{filter_method = "mixture"}: \code{"optimal"} (default, precomputed
#'   tables — deterministic), \code{"ksc"} (published KSC table; Gaussian,
#'   \code{K = 7} only) or \code{"em"} (legacy sampled EM fit; any \code{K}).
#'   See \code{\link{filter_svp}} for details.
#' @param K Integer. Number of mixture components for GMKF. Default 7.
#'   Under \code{mixture = "optimal"} only 7, 10, 12 or 16 are available.
#' @param M Integer. Number of particles for BPF. Default 1000.
#' @param seed Integer. Random seed for BPF and the \code{mixture = "em"}
#'   fit. Default 42.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#'
#' @details
#' The forecasts follow Ahsan, Dufour and Rodriguez-Rondon (2026), Section
#' 3.3.6. Writing \eqn{\hat{\xi}_{T|T}} and \eqn{P_{T|T}} for the filtered state
#' and covariance at \eqn{T} (from \code{\link{filter_svp}}),
#' \deqn{\hat{\xi}_{T+h|T} = F^h \hat{\xi}_{T|T} + F^{h-1} r\, \sigma_v \delta_p
#'   \hat{z}_T, \qquad P_{T+h|T} = F P_{T+h-1|T} F' + \sigma_v^2 r r' \ (h \ge 2),}
#' where \eqn{\hat{z}_T = g(\hat{u}_T)} is the Bayes-optimal leverage proxy of
#' the last filtered residual and \eqn{P_{T+1|T}} is the one-step prediction
#' covariance of eqs. (48a)-(48b), \eqn{A_{T+1} P_{T|T} A_{T+1}' + \sigma_v^2
#' [(1-\delta_p^2) + \delta_p^2 V_T] r r'} with \eqn{A_{T+1}} built from
#' \eqn{\hat{u}_T} (see \code{\link{filter_svp}} for \eqn{g}, \eqn{A_t} and
#' \eqn{V_t}). Beyond one step no further leverage shock is observed, so the
#' full marginal state-innovation variance \eqn{\sigma_v^2} applies. Without
#' leverage the recursion is the standard \eqn{P_{T+h|T} = F P_{T+h-1|T} F' +
#' \sigma_v^2 r r'} at every horizon. For the particle filter
#' \eqn{P_{T|T}} is the diagonal approximation returned by
#' \code{\link{filter_svp}}.
#'
#' @return An object of class \code{"svp_forecast"}, a list containing:
#' \describe{
#'   \item{w_forecasted}{Primary forecast (scale determined by \code{output}).}
#'   \item{log_var_forecast}{Log-volatility forecasts \eqn{w_{T+h|T}}.}
#'   \item{var_forecast}{Conditional variance forecasts \eqn{\sigma^2_{T+h|T}}.}
#'   \item{vol_forecast}{Conditional volatility forecasts \eqn{\sigma_{T+h|T}}.}
#'   \item{P_forecast}{Forecast MSE \eqn{P_{T+h|T}} for each horizon.}
#'   \item{w_estimated}{Filtered log-volatility.}
#'   \item{w_smoothed}{Smoothed log-volatility.}
#'   \item{zt}{Filtered standardized residuals.}
#'   \item{zt_smoothed}{Smoothed standardized residuals.}
#'   \item{ys}{Demeaned log-squared returns.}
#'   \item{mdl}{The estimated model object.}
#'   \item{H}{The forecast horizon.}
#'   \item{output}{The chosen output scale.}
#'   \item{filter_output}{The \code{"svp_filter"} object from filtering.}
#' }
#'
#' @examples
#' \donttest{
#' sim <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                leverage = TRUE, rho = -0.3)
#' fit <- svp(sim$y, p = 1, leverage = TRUE)
#' fc <- forecast_svp(fit, H = 10)
#' plot(fc)
#' }
#'
#' @export
forecast_svp <- function(object, H = 1,
                         output = c("log-variance", "variance", "volatility"),
                         filter_method = "corrected",
                         mixture = c("optimal", "ksc", "em"),
                         K = 7, M = 1000, seed = 42, del = 1e-10) {
  if (!inherits(object, c("svp", "svp_t", "svp_ged"))) {
    stop("object must be an svp/svp_t/svp_ged model from svp(). ",
         "Usage: fit <- svp(y, ...); fc <- forecast_svp(fit, H = 10)")
  }
  if (!is.numeric(H) || length(H) != 1L || H < 1L)
    stop("'H' must be a positive integer (>= 1).")
  H <- as.integer(H)
  if (!is.numeric(K) || length(K) != 1L || K < 1L)
    stop("'K' must be a positive integer (>= 1).")
  if (!is.numeric(M) || length(M) != 1L || M < 1L)
    stop("'M' must be a positive integer (>= 1).")
  if (!is.numeric(del) || length(del) != 1L || del <= 0)
    stop("'del' must be a positive number.")
  output <- match.arg(output)
  filter_method <- match.arg(filter_method, c("corrected", "mixture", "particle"))
  mixture <- match.arg(mixture)

  mdl <- object
  y_vec <- as.numeric(mdl$y)
  phi <- mdl$phi
  p_len <- length(phi)
  delta_p <- if (is.null(mdl$rho) || is.na(mdl$rho)) 0 else mdl$rho
  sigma_y <- mdl$sigy
  sigma_v <- mdl$sigv

  # Filter
  filt <- filter_svp(mdl, method = filter_method,
                     mixture = mixture, K = K, M = M, seed = seed, del = del)
  Tsize <- length(filt$w_filtered)

  # Final filtered state and standardized residual u_hat_T
  xi_T <- filt$xi_filtered[, Tsize]
  uT <- filt$zt[Tsize]
  # Leverage proxy g(u_hat_T) = E[z_T | u_T] with its derivative and the
  # mixing-variable variance Var(z_T | u_T) (R mirror of the C++ helper)
  lp <- .leverage_proxy_R(uT, mdl)
  zT <- lp$g

  # Build companion matrix
  Fmat <- companionMat(t(phi), p_len, 1)
  h_vec <- c(1, rep(0, p_len - 1))
  r_vec <- h_vec
  rr <- r_vec %*% t(r_vec)

  # Prediction covariance, Ahsan, Dufour and Rodriguez-Rondon (2026) Section 3.3.6.
  # h = 1: the one-step prediction of eqs. (48a)-(48b) from (xi_T, P_T, u_hat_T):
  #          P_{T+1|T} = A_1 P_T A_1' + sigma_v^2 [(1 - delta^2) + delta^2 V_T] r r',
  #          A_1 = F - (1/2) sigma_v delta g'(u_hat_T) u_hat_T r h'.
  # h >= 2: no further leverage shock is observed, so the marginal
  #          state-innovation variance applies: P_h = F P_{h-1} F' + sigma_v^2 r r'.
  # delta = 0 gives A_1 = F and sigma_v^2 at every horizon (no-leverage case).
  if (delta_p != 0) {
    A_1 <- Fmat
    A_1[1, 1] <- A_1[1, 1] - 0.5 * sigma_v * delta_p * lp$dg * uT
    q_1 <- sigma_v^2 * ((1 - delta_p^2) + delta_p^2 * lp$V)
  } else {
    A_1 <- Fmat
    q_1 <- sigma_v^2
  }
  Q_hge2 <- sigma_v^2 * rr

  # E[u^2] for variance forecast
  Eu2 <- 1.0
  if (inherits(mdl, "svp_t") && !is.null(mdl$v) && is.finite(mdl$v) && mdl$v > 2) {
    Eu2 <- mdl$v / (mdl$v - 2)
  }

  # h-step forecasts
  tilde_delta_p <- delta_p * sigma_v
  eta_T <- c(zT, rep(0, p_len - 1))

  # Get final filtered covariance for MSE recursion
  # Use full p x p matrix from filter output (exact for p>=2)
  P_T <- filt$P_filt_T

  log_var_fc <- numeric(H)
  P_fc <- numeric(H)
  var_fc <- numeric(H)
  vol_fc <- numeric(H)

  xi_h <- xi_T
  P_h <- P_T

  for (h in 1:H) {
    # State prediction (leverage shift enters at h = 1 only)
    if (h == 1) {
      xi_h <- as.numeric(Fmat %*% xi_T) + tilde_delta_p * eta_T
      P_h <- A_1 %*% P_h %*% t(A_1) + q_1 * rr
    } else {
      xi_h <- as.numeric(Fmat %*% xi_h)
      P_h <- Fmat %*% P_h %*% t(Fmat) + Q_hge2
    }

    w_h <- sum(h_vec * xi_h)
    P_h_11 <- as.numeric(t(h_vec) %*% P_h %*% h_vec)

    log_var_fc[h] <- w_h
    P_fc[h] <- P_h_11
    # Variance forecast: sigma_y^2 * E[u^2] * E[exp(w)] where w ~ N(w_h, P_h_11)
    var_fc[h] <- sigma_y^2 * Eu2 * exp(w_h + P_h_11 / 2)
    vol_fc[h] <- sqrt(var_fc[h])
  }

  # Primary output
  w_forecasted <- switch(output,
    "log-variance" = log_var_fc,
    "variance"     = var_fc,
    "volatility"   = vol_fc
  )

  # Build output (backward compatible fields)
  ys <- log(y_vec^2 + del) - mdl$mu
  out <- list(
    w_forecasted    = w_forecasted,
    log_var_forecast = log_var_fc,
    var_forecast    = var_fc,
    vol_forecast    = vol_fc,
    P_forecast      = P_fc,
    w_estimated     = matrix(filt$w_filtered, ncol = 1),
    w_smoothed      = matrix(filt$w_smoothed, ncol = 1),
    zt              = matrix(filt$zt, ncol = 1),
    zt_smoothed     = matrix(filt$zt_smoothed, ncol = 1),
    ys              = as.numeric(ys),
    mdl             = mdl,
    H               = H,
    output          = output,
    filter_output   = filt,
    call            = match.call()
  )
  class(out) <- "svp_forecast"
  return(out)
}


# Internal: R-side mirror of the C++ helpers leverage_proxy (kalman_filter.cpp)
# and ged_copula_cpp (utils_cpp.cpp).  For the standardized residual u returns
#   g  = E[z | u]        the Bayes-optimal leverage proxy (SVHT Section 3.3.2),
#   dg = g'(u)           for the Jacobian A_t of eq. (48b),
#   V  = Var(z | u)      for the state-noise variance of eq. (48a).
# Gaussian: g = u, dg = 1, V = 0.  GED: g = Phi^{-1}(F_GED(u; nu)) (eq. 49),
# computed tail-stably from the log upper-tail probability exactly as in C++
# (R's pgamma/qnorm are the same routines), dg = f_GED(u)/phi(g), V = 0.
# Student-t: c = Gamma((nu+2)/2)/Gamma((nu+1)/2), g = u sqrt(2/(nu+u^2)) c,
# dg = sqrt(2) nu c (nu+u^2)^{-3/2}, V = u^2 [(nu+1) - 2 c^2]/(nu+u^2).
#' @keywords internal
.leverage_proxy_R <- function(u, fit) {
  if (inherits(fit, "svp")) return(list(g = u, dg = 1, V = 0))
  nu <- as.numeric(fit$v)
  if (inherits(fit, "svp_ged")) {
    a <- exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu)))
    x <- (abs(u) / a)^nu
    log_q <- log(0.5) + pgamma(x, shape = 1 / nu, lower.tail = FALSE, log.p = TRUE)
    g_abs <- qnorm(log_q, lower.tail = FALSE, log.p = TRUE)
    g <- if (u < 0) -g_abs else g_abs
    log_f <- log(nu) - log(2 * a) - lgamma(1 / nu) - x
    return(list(g = g, dg = exp(log_f - dnorm(g, log = TRUE)), V = 0))
  }
  # Student-t (svp_t)
  cc <- exp(lgamma((nu + 2) / 2) - lgamma((nu + 1) / 2))
  s <- nu + u^2
  list(g  = u * sqrt(2 / s) * cc,
       dg = sqrt(2) * nu * cc / (s * sqrt(s)),
       V  = u^2 * ((nu + 1) / s - 2 * cc^2 / s))
}
