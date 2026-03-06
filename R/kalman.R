#' Kalman Filter and Smoother for SV(p) Models
#'
#' Applies Kalman filtering and Rauch-Tung-Striebel smoothing to estimate the
#' latent log-volatility process from an estimated SV(p) model.
#'
#' @param y Numeric vector. Observed returns.
#' @param model An \code{"svp"} object from \code{\link{svp}}.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#'
#' @return A list with:
#' \describe{
#'   \item{w_estimated}{Filtered log-volatility (Kalman filter output).}
#'   \item{w_smoothed}{Smoothed log-volatility (Kalman smoother output).}
#'   \item{zt}{Filtered standardized residuals.}
#'   \item{zt_smoothed}{Smoothed standardized residuals.}
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
  if (!inherits(model, "svp")) stop("model must be of class 'svp'.")
  y <- as.matrix(as.numeric(y))
  mu <- model$mu
  phi <- model$phi
  p <- length(phi)
  delta_p <- if (is.null(model$rho) || is.na(model$rho)) 0 else model$rho
  sigma_y <- model$sigy
  sigma_v <- model$sigv

  ys <- log(y^2 + del) - mu
  Tsize <- length(ys)

  Fmat <- companionMat(t(phi), p, 1)
  Hprm <- t(as.matrix(c(1, rep(0, p - 1))))
  etat <- as.matrix(c(1, rep(0, p - 1)))
  Q <- matrix(0, p, p)
  Q[1, 1] <- 1

  zt <- matrix(0, Tsize, 1)
  w_estimated <- matrix(0, Tsize, 1)
  xi_pred <- matrix(0, p, Tsize)
  xi_estimated <- matrix(0, p, Tsize)
  P_pred <- array(0, c(p, p, Tsize))
  P_updated <- array(0, c(p, p, Tsize))

  # Initialize
  xi_pred[, 1] <- rep(0, p)
  P_pred[, , 1] <- diag(p) * (sigma_v^2)

  # Forward pass (filtering)
  for (t in 1:Tsize) {
    # Update
    K <- P_pred[, , t] %*% t(Hprm) %*%
      (1 / (Hprm %*% P_pred[, , t] %*% t(Hprm) + (pi^2) / 2))
    xi_estimated[, t] <- xi_pred[, t, drop = FALSE] +
      K %*% (ys[t] - Hprm %*% xi_pred[, t, drop = FALSE])
    P_updated[, , t] <- P_pred[, , t] - K %*% Hprm %*% P_pred[, , t]
    w_estimated[t, ] <- Hprm %*% xi_estimated[, t, drop = FALSE]
    zt[t, ] <- (y[t, ] / sigma_y) * exp((-1 / 2) * w_estimated[t, ])
    # Predict
    if (t < Tsize) {
      xi_pred[, t + 1] <- Fmat %*% xi_estimated[, t, drop = FALSE] +
        (delta_p * sigma_v) * (etat * zt[t, ])
      P_pred[, , t + 1] <- Fmat %*% P_updated[, , t] %*% t(Fmat) +
        (sigma_v^2) * (delta_p^2) * Q + (sigma_v^2) * (1 - (delta_p^2)) * Q
    }
  }

  # Backward pass (smoothing)
  w_smoothed <- matrix(0, Tsize, 1)
  xi_smoothed <- matrix(0, p, Tsize)
  xi_smoothed[, Tsize] <- xi_estimated[, Tsize]
  w_smoothed[Tsize, ] <- Hprm %*% xi_estimated[, Tsize]
  P_smoothed <- array(0, c(p, p, Tsize))
  P_smoothed[, , Tsize] <- P_updated[, , Tsize]

  for (t in (Tsize - 1):1) {
    A <- P_updated[, , t] %*% t(Fmat) %*% solve(P_pred[, , t + 1])
    xi_smoothed[, t] <- xi_estimated[, t] +
      A %*% (xi_smoothed[, t + 1] - xi_pred[, t + 1])
    w_smoothed[t, ] <- Hprm %*% xi_smoothed[, t]
    P_smoothed[, , t] <- P_updated[, , t] +
      A %*% (P_smoothed[, , t + 1] - P_pred[, , t + 1]) %*% t(A)
  }

  zt_smoothed <- (y / sigma_y) * exp((-1 / 2) * w_smoothed)

  list(w_estimated = w_estimated, w_smoothed = w_smoothed,
       zt = zt, zt_smoothed = zt_smoothed)
}
