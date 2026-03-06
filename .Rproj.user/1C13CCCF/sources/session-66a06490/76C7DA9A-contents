# =========================================================================== #
# S3 methods for wARMASVp classes
# =========================================================================== #

# ----------- svp class methods ----------- #

#' @export
print.svp <- function(x, ...) {
  p <- x$p
  mdl_type <- if (isTRUE(x$leverage) && !is.na(x$rho)) "SVL" else "SV"
  cat(sprintf("%s(%d) model (W-ARMA-SV estimation)\n", mdl_type, p))
  cat(sprintf("  phi: %s\n", paste(round(x$phi, 4), collapse = ", ")))
  cat(sprintf("  sigy: %.4f, sigv: %.4f\n", x$sigy, x$sigv))
  if (!is.na(x$rho)) {
    cat(sprintf("  rho (leverage): %.4f\n", x$rho))
  }
  cat(sprintf("  T = %d, J = %d\n", length(x$y), x$J))
  invisible(x)
}

#' @export
summary.svp <- function(object, ...) {
  p <- object$p
  mdl_type <- if (isTRUE(object$leverage) && !is.na(object$rho)) "SVL" else "SV"
  cat(sprintf("\n%s(%d) Model - W-ARMA-SV Estimation\n", mdl_type, p))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Sample size: %d\n", length(object$y)))
  cat(sprintf("Winsorizing parameter J: %d\n", object$J))
  if (!is.na(object$rho)) {
    cat(sprintf("Leverage correlation type: %s\n", object$rho_type))
  }
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Parameter estimates:\n\n")
  # Build parameter table
  pnames <- c(paste0("phi_", 1:p), "sigma_y", "sigma_v")
  pvals <- c(object$phi, object$sigy, object$sigv)
  if (!is.na(object$rho)) {
    pnames <- c(pnames, "rho")
    pvals <- c(pvals, object$rho)
  }
  df <- data.frame(Parameter = pnames, Estimate = round(pvals, 6))
  print(df, row.names = FALSE)
  cat("\n")
  if (!is.na(object$gammatilde)) {
    cat(sprintf("gamma_tilde: %.6f\n", object$gammatilde))
  }
  if (isTRUE(object$nonstationary_ind)) {
    cat("Note: Stationarity was enforced (roots projected inside unit circle).\n")
  }
  cat("\n")
  invisible(object)
}

#' @export
coef.svp <- function(object, ...) {
  p <- object$p
  nms <- c(paste0("phi_", 1:p), "sigma_y", "sigma_v")
  vals <- c(object$phi, object$sigy, object$sigv)
  if (!is.na(object$rho)) {
    nms <- c(nms, "rho")
    vals <- c(vals, object$rho)
  }
  names(vals) <- nms
  return(vals)
}

# ----------- svp_t class methods ----------- #

#' @export
print.svp_t <- function(x, ...) {
  cat("SV(1) model with Student-t errors (W-ARMA-SV estimation)\n")
  cat(sprintf("  phi: %.4f\n", x$phi))
  cat(sprintf("  sigy: %.4f, sigv: %.4f\n", x$sigy, x$sigv))
  cat(sprintf("  nu (d.f.): %.4f\n", x$v))
  cat(sprintf("  T = %d, J = %d\n", length(x$y), x$J))
  invisible(x)
}

#' @export
summary.svp_t <- function(object, ...) {
  cat("\nSV(1) Model with Student-t Errors\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Sample size: %d\n", length(object$y)))
  cat(sprintf("Winsorizing parameter J: %d\n", object$J))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Parameter estimates:\n\n")
  pnames <- c("phi", "sigma_y", "sigma_v", "nu")
  pvals <- c(object$phi, object$sigy, object$sigv, object$v)
  df <- data.frame(Parameter = pnames, Estimate = round(pvals, 6))
  print(df, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' @export
coef.svp_t <- function(object, ...) {
  vals <- c(object$phi, object$sigy, object$sigv, object$v)
  names(vals) <- c("phi", "sigma_y", "sigma_v", "nu")
  return(vals)
}

# ----------- svp_ged class methods ----------- #

#' @export
print.svp_ged <- function(x, ...) {
  cat("SV(1) model with GED errors (W-ARMA-SV estimation)\n")
  cat(sprintf("  phi: %.4f\n", x$phi))
  cat(sprintf("  sigy: %.4f, sigv: %.4f\n", x$sigy, x$sigv))
  cat(sprintf("  nu (shape): %.4f\n", x$v))
  cat(sprintf("  T = %d, J = %d\n", length(x$y), x$J))
  invisible(x)
}

#' @export
summary.svp_ged <- function(object, ...) {
  cat("\nSV(1) Model with GED Errors\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Sample size: %d\n", length(object$y)))
  cat(sprintf("Winsorizing parameter J: %d\n", object$J))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Parameter estimates:\n\n")
  pnames <- c("phi", "sigma_y", "sigma_v", "nu")
  pvals <- c(object$phi, object$sigy, object$sigv, object$v)
  df <- data.frame(Parameter = pnames, Estimate = round(pvals, 6))
  print(df, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' @export
coef.svp_ged <- function(object, ...) {
  vals <- c(object$phi, object$sigy, object$sigv, object$v)
  names(vals) <- c("phi", "sigma_y", "sigma_v", "nu")
  return(vals)
}

# ----------- svp_forecast class methods ----------- #

#' @export
print.svp_forecast <- function(x, ...) {
  H <- x$H
  mdl <- x$mdl
  p <- mdl$p
  mdl_type <- if (isTRUE(mdl$leverage) && !is.na(mdl$rho)) "SVL" else "SV"
  cat(sprintf("%s(%d) Forecast (H = %d)\n", mdl_type, p, H))
  cat(sprintf("  Forecasted log-volatility:\n"))
  for (h in seq_len(min(H, 10))) {
    cat(sprintf("    h=%d: %.4f\n", h, x$w_forecasted[h]))
  }
  if (H > 10) cat(sprintf("    ... (%d more horizons)\n", H - 10))
  invisible(x)
}

#' @export
plot.svp_forecast <- function(x, ...) {
  Tsize <- length(x$w_smoothed)
  H <- x$H
  all_w <- c(x$w_smoothed, x$w_forecasted)
  time_idx <- seq_len(Tsize + H)
  plot(time_idx, all_w, type = "n",
       xlab = "Time", ylab = "Log-volatility",
       main = "Smoothed + Forecasted Log-Volatility", ...)
  graphics::lines(1:Tsize, x$w_smoothed, col = "black")
  graphics::lines((Tsize + 1):(Tsize + H), x$w_forecasted, col = "red", lwd = 2)
  graphics::abline(v = Tsize + 0.5, lty = 2, col = "gray50")
  graphics::legend("topright", legend = c("Smoothed", "Forecast"),
                   col = c("black", "red"), lty = 1, lwd = c(1, 2))
  invisible(x)
}

# ----------- svp_test class methods ----------- #

#' @export
print.svp_test <- function(x, ...) {
  cat(sprintf("%s Test\n", x$test_type))
  cat(paste(rep("-", 40), collapse = ""), "\n")
  for (i in seq_along(x$null_param)) {
    cat(sprintf("  H0: %s = %s\n", x$null_param[i], x$null_value[i]))
  }
  cat(sprintf("  Test statistic: %.4f\n", x$s0))
  cat(sprintf("  p-value: %.4f\n", x$pval))
  cat(sprintf("  MC replications: %d\n", length(x$sN)))
  cat("\n")
  invisible(x)
}
