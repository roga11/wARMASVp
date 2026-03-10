#' Multi-Step Ahead Volatility Forecast
#'
#' Estimates an SV(p) model, applies Kalman filtering/smoothing, and produces
#' multi-step ahead forecasts of log-volatility.
#'
#' @param y Numeric vector. Observed returns.
#' @param p Integer. Order of the volatility process. Default 1.
#' @param J Integer. Winsorizing parameter. Default 10.
#' @param H Integer. Maximum forecast horizon. Default 1.
#' @param leverage Logical. Estimate leverage. Default \code{TRUE}.
#' @param errorType Character. Innovation distribution: \code{"Gaussian"},
#'   \code{"Student-t"}, or \code{"GED"}. Default \code{"Gaussian"}.
#' @param rho_type Character. Correlation type for leverage. Default \code{"pearson"}.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#' @param trunc_lev Logical. Truncate leverage correlation estimate to
#'   \code{[-1,1]}. Default \code{TRUE}.
#' @param wDecay Logical. Use decaying weights. Default \code{FALSE}.
#' @param logNu Logical. Solve for \eqn{\nu} in log-space for numerical
#'   stability (Student-t only). Default \code{FALSE}.
#'
#' @return An object of class \code{"svp_forecast"}, a list containing:
#' \describe{
#'   \item{w_estimated}{Filtered log-volatility.}
#'   \item{w_smoothed}{Smoothed log-volatility.}
#'   \item{w_forecasted}{Forecasted log-volatility for horizons 1 to H.}
#'   \item{zt}{Filtered standardized residuals.}
#'   \item{zt_smoothed}{Smoothed standardized residuals.}
#'   \item{ys}{Demeaned log-squared returns.}
#'   \item{mdl}{The estimated \code{"svp"} model object.}
#'   \item{H}{The forecast horizon.}
#' }
#'
#' @note The Kalman filter update step uses a Gaussian approximation for the
#'   measurement equation regardless of \code{errorType}. When
#'   \code{errorType} is \code{"Student-t"} or \code{"GED"}, the
#'   distribution-specific measurement noise variance \eqn{\sigma_\varepsilon^2}
#'   is used, but the filter update remains linear-Gaussian. This is standard
#'   practice for non-Gaussian state-space models.
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, leverage = TRUE, rho = -0.3)$y
#' fc <- forecast_svp(y, p = 1, H = 10)
#' plot(fc$w_forecasted, type = "h", main = "Forecasted log-volatility")
#' }
#'
#' @importFrom expm %^%
#' @export
forecast_svp <- function(y, p = 1, J = 10, H = 1, leverage = TRUE,
                         errorType = "Gaussian", rho_type = "pearson",
                         del = 1e-10, trunc_lev = TRUE, wDecay = FALSE,
                         logNu = FALSE) {
  y <- as.matrix(as.numeric(y))
  mdl <- svp(as.numeric(y), p, J, leverage = leverage, errorType = errorType,
             rho_type = rho_type, del = del, trunc_lev = trunc_lev,
             wDecay = wDecay, logNu = logNu)

  mu <- mdl$mu
  phi <- mdl$phi
  p_len <- length(phi)
  delta_p <- if (is.null(mdl$rho) || is.na(mdl$rho)) 0 else mdl$rho
  sigma_y <- mdl$sigy
  sigma_v <- mdl$sigv

  ys <- log(y^2 + del) - mu
  kalout <- kalman_filter(as.numeric(y), mdl, del)
  w_estimated <- kalout$w_smoothed
  zt <- kalout$zt
  Tsize <- length(ys)

  tilde_delta_p <- delta_p * sigma_v
  Fmat <- companionMat(t(phi), p_len, 1)
  Hprm <- t(as.matrix(c(1, rep(0, p_len - 1))))
  ztildevec <- as.matrix(c(1, rep(0, p_len - 1)))

  zT <- zt[Tsize, ]
  xi_T <- as.matrix(w_estimated[Tsize:(Tsize - p_len + 1)])
  eta_T <- as.matrix(c(zT, rep(0, p_len - 1)))

  w_forecasted <- numeric(H)
  for (h in 1:H) {
    Fh <- expm::`%^%`(Fmat, h)
    Fh1 <- if (h == 1) diag(p_len) else expm::`%^%`(Fmat, h - 1)
    w_forecasted[h] <- Hprm %*% Fh %*% xi_T +
      tilde_delta_p * Hprm %*% Fh1 %*% eta_T
  }

  out <- kalout
  out$w_forecasted <- w_forecasted - del
  out$ys <- as.numeric(ys)
  out$mdl <- mdl
  out$H <- H
  out$call <- match.call()
  class(out) <- "svp_forecast"
  return(out)
}
