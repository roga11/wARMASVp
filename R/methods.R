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
  p <- length(x$phi)
  mdl_type <- if (isTRUE(x$leverage) && !is.na(x$rho)) "SVL" else "SV"
  cat(sprintf("%s(%d) model with Student-t errors (W-ARMA-SV estimation)\n", mdl_type, p))
  cat(sprintf("  phi: %s\n", paste(round(x$phi, 4), collapse = ", ")))
  cat(sprintf("  sigy: %.4f, sigv: %.4f\n", x$sigy, x$sigv))
  cat(sprintf("  nu (d.f.): %.4f\n", x$v))
  if (!is.na(x$rho)) {
    cat(sprintf("  rho (leverage): %.4f\n", x$rho))
  }
  cat(sprintf("  T = %d, J = %d\n", length(x$y), x$J))
  invisible(x)
}

#' @export
summary.svp_t <- function(object, ...) {
  p <- length(object$phi)
  mdl_type <- if (isTRUE(object$leverage) && !is.na(object$rho)) "SVL" else "SV"
  cat(sprintf("\n%s(%d) Model with Student-t Errors\n", mdl_type, p))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Sample size: %d\n", length(object$y)))
  cat(sprintf("Winsorizing parameter J: %d\n", object$J))
  if (!is.na(object$rho)) {
    cat(sprintf("Leverage correlation type: %s\n", object$rho_type))
  }
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Parameter estimates:\n\n")
  pnames <- c(paste0("phi_", seq_len(p)), "sigma_y", "sigma_v", "nu")
  pvals <- c(object$phi, object$sigy, object$sigv, object$v)
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
  cat("\n")
  invisible(object)
}

#' @export
coef.svp_t <- function(object, ...) {
  p <- length(object$phi)
  nms <- c(paste0("phi_", seq_len(p)), "sigma_y", "sigma_v", "nu")
  vals <- c(object$phi, object$sigy, object$sigv, object$v)
  if (!is.na(object$rho)) {
    nms <- c(nms, "rho")
    vals <- c(vals, object$rho)
  }
  names(vals) <- nms
  return(vals)
}

# ----------- svp_ged class methods ----------- #

#' @export
print.svp_ged <- function(x, ...) {
  p <- length(x$phi)
  mdl_type <- if (isTRUE(x$leverage) && !is.na(x$rho)) "SVL" else "SV"
  cat(sprintf("%s(%d) model with GED errors (W-ARMA-SV estimation)\n", mdl_type, p))
  cat(sprintf("  phi: %s\n", paste(round(x$phi, 4), collapse = ", ")))
  cat(sprintf("  sigy: %.4f, sigv: %.4f\n", x$sigy, x$sigv))
  cat(sprintf("  nu (shape): %.4f\n", x$v))
  if (!is.na(x$rho)) {
    cat(sprintf("  rho (leverage): %.4f\n", x$rho))
  }
  cat(sprintf("  T = %d, J = %d\n", length(x$y), x$J))
  invisible(x)
}

#' @export
summary.svp_ged <- function(object, ...) {
  p <- length(object$phi)
  mdl_type <- if (isTRUE(object$leverage) && !is.na(object$rho)) "SVL" else "SV"
  cat(sprintf("\n%s(%d) Model with GED Errors\n", mdl_type, p))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Sample size: %d\n", length(object$y)))
  cat(sprintf("Winsorizing parameter J: %d\n", object$J))
  if (!is.na(object$rho)) {
    cat(sprintf("Leverage correlation type: %s\n", object$rho_type))
  }
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Parameter estimates:\n\n")
  pnames <- c(paste0("phi_", seq_len(p)), "sigma_y", "sigma_v", "nu")
  pvals <- c(object$phi, object$sigy, object$sigv, object$v)
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
  cat("\n")
  invisible(object)
}

#' @export
coef.svp_ged <- function(object, ...) {
  p <- length(object$phi)
  nms <- c(paste0("phi_", seq_len(p)), "sigma_y", "sigma_v", "nu")
  vals <- c(object$phi, object$sigy, object$sigv, object$v)
  if (!is.na(object$rho)) {
    nms <- c(nms, "rho")
    vals <- c(vals, object$rho)
  }
  names(vals) <- nms
  return(vals)
}

# ----------- svp_filter class methods ----------- #

#' @export
print.svp_filter <- function(x, ...) {
  T_obs <- length(x$w_filtered)
  mdl <- x$model
  p <- length(mdl$phi)
  mdl_type <- if (isTRUE(mdl$leverage) && !is.na(mdl$rho)) "SVL" else "SV"
  dist_label <- if (inherits(mdl, "svp_t")) "Student-t"
                else if (inherits(mdl, "svp_ged")) "GED"
                else "Gaussian"
  cat(sprintf("%s(%d) Filter [%s, %s]\n", mdl_type, p, x$method, dist_label))
  cat(sprintf("  T = %d, log-likelihood = %.2f\n", T_obs, x$loglik))
  cat(sprintf("  Mean filtered MSE: %.4f\n", mean(x$P_filtered)))
  invisible(x)
}

#' @export
summary.svp_filter <- function(object, ...) {
  print.svp_filter(object, ...)
  cat(sprintf("  w_filtered range: [%.4f, %.4f]\n",
              min(object$w_filtered), max(object$w_filtered)))
  cat(sprintf("  w_smoothed range: [%.4f, %.4f]\n",
              min(object$w_smoothed), max(object$w_smoothed)))
  invisible(object)
}

#' @export
plot.svp_filter <- function(x, ...) {
  T_obs <- length(x$w_filtered)
  time_idx <- seq_len(T_obs)
  ylim <- range(c(x$w_filtered, x$w_smoothed))
  plot(time_idx, x$w_filtered, type = "l", col = "steelblue",
       xlab = "Time", ylab = "Log-volatility",
       main = sprintf("Filtered & Smoothed Log-Volatility [%s]", x$method),
       ylim = ylim, ...)
  graphics::lines(time_idx, x$w_smoothed, col = "darkred", lwd = 1.5)
  graphics::legend("topright", legend = c("Filtered", "Smoothed"),
                   col = c("steelblue", "darkred"), lty = 1, lwd = c(1, 1.5))
  invisible(x)
}

# ----------- svp_forecast class methods ----------- #

#' @export
print.svp_forecast <- function(x, ...) {
  H <- x$H
  mdl <- x$mdl
  p <- mdl$p
  mdl_type <- if (isTRUE(mdl$leverage) && !is.na(mdl$rho)) "SVL" else "SV"
  out_label <- if (!is.null(x$output)) x$output else "log-variance"
  cat(sprintf("%s(%d) Forecast (H = %d, output = %s)\n", mdl_type, p, H, out_label))
  cat(sprintf("  Forecasted values:\n"))
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
  out_label <- if (!is.null(x$output)) x$output else "log-variance"
  # For plot, always show log-variance (smoothed + forecast)
  fc_vals <- if (!is.null(x$log_var_forecast)) x$log_var_forecast else x$w_forecasted
  all_w <- c(x$w_smoothed, fc_vals)
  time_idx <- seq_len(Tsize + H)
  plot(time_idx, all_w, type = "n",
       xlab = "Time", ylab = "Log-volatility",
       main = "Smoothed + Forecasted Log-Volatility", ...)
  graphics::lines(1:Tsize, x$w_smoothed, col = "black")
  graphics::lines((Tsize + 1):(Tsize + H), fc_vals, col = "red", lwd = 2)
  graphics::abline(v = Tsize + 0.5, lty = 2, col = "gray50")
  graphics::legend("topright", legend = c("Smoothed", "Forecast"),
                   col = c("black", "red"), lty = 1, lwd = c(1, 2))
  invisible(x)
}

# ----------- svp_test class methods ----------- #

#' @export
print.svp_test <- function(x, ...) {
  dir_label <- if (!is.null(x$direction) && x$direction != "two-sided") {
    sprintf(" [%s]", x$direction)
  } else ""
  cat(sprintf("%s Test%s\n", x$test_type, dir_label))
  cat(paste(rep("-", 40), collapse = ""), "\n")
  for (i in seq_along(x$null_param)) {
    cat(sprintf("  H0: %s = %s\n", x$null_param[i], x$null_value[i]))
  }
  cat(sprintf("  Test statistic (LR): %.4f\n", x$s0))
  if (!is.null(x$S_T)) {
    cat(sprintf("  Signed root (S_T): %.4f\n", x$S_T))
  }
  cat(sprintf("  p-value: %.4f\n", x$pval))
  cat(sprintf("  MC replications: %d\n", length(x$sN)))
  cat("\n")
  invisible(x)
}
