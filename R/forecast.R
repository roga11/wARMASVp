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
#' @param K Integer. Number of mixture components for GMKF. Default 7.
#' @param M Integer. Number of particles for BPF. Default 1000.
#' @param seed Integer. Random seed for BPF. Default 42.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
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

  mdl <- object
  y_vec <- as.numeric(mdl$y)
  phi <- mdl$phi
  p_len <- length(phi)
  delta_p <- if (is.null(mdl$rho) || is.na(mdl$rho)) 0 else mdl$rho
  sigma_y <- mdl$sigy
  sigma_v <- mdl$sigv

  # Filter
  filt <- filter_svp(mdl, method = filter_method, K = K, M = M, seed = seed,
                     del = del)
  Tsize <- length(filt$w_filtered)

  # Extract final state and residual
  xi_T <- filt$xi_filtered[, Tsize]
  zT <- filt$zt[Tsize]

  # Build companion matrix
  Fmat <- companionMat(t(phi), p_len, 1)
  h_vec <- c(1, rep(0, p_len - 1))
  r_vec <- h_vec

  # State noise covariance (with leverage component for Student-t)
  # For Gaussian/GED: var_zt=1, so Q = sigma_v^2 (unchanged)
  # For Student-t leverage: var_zt=nu/(nu-2), Q > sigma_v^2
  var_zt_fc <- 1.0
  if (inherits(mdl, "svp_t") && delta_p != 0 &&
      !is.null(mdl$v) && is.finite(mdl$v) && mdl$v > 2) {
    var_zt_fc <- mdl$v / (mdl$v - 2)
  }
  Q_scalar <- if (delta_p == 0) sigma_v^2
              else sigma_v^2 * (1 - delta_p^2 + delta_p^2 * var_zt_fc)
  Q <- Q_scalar * (r_vec %*% t(r_vec))

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
    # State prediction (no leverage feedback beyond h=1)
    if (h == 1) {
      xi_h <- as.numeric(Fmat %*% xi_T) + tilde_delta_p * eta_T
    } else {
      xi_h <- as.numeric(Fmat %*% xi_h)
    }
    P_h <- Fmat %*% P_h %*% t(Fmat) + Q

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
