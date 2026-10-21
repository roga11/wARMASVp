# Tests for the leverage prediction covariance of the CKF / GMKF and of the
# h = 1 forecast (Ahsan, Dufour and Rodriguez-Rondon 2026, eqs. 48a-48b and
# Section 3.3.6), introduced in 0.2.0.9006, which also removed the `proxy`
# argument (the Bayes-optimal proxy is the only one). Earlier filter fixes
# kept here: BPF posterior lambda for Student-t leverage (2026-05-09, Fix 3).

# R reference for the CKF under leverage at p = 1: SVHT (48a)-(48b) with the
# per-family g, g', V written out. Independent of the package filter code.
.ckf_ref <- function(fit, y, del = 1e-10) {
  pr <- wARMASVp:::.get_filter_params(fit)
  nu <- pr$nu; sv <- pr$sigma_v; d <- pr$delta_p; phi <- as.numeric(pr$phi)
  gdv <- switch(pr$dist_name,
    gaussian  = function(u) c(u, 1, 0),
    student_t = function(u) {
      cc <- exp(lgamma((nu + 2) / 2) - lgamma((nu + 1) / 2)); s <- nu + u^2
      c(u * sqrt(2 / s) * cc, sqrt(2) * nu * cc / s^1.5,
        u^2 * ((nu + 1) / s - 2 * cc^2 / s)) },
    ged       = function(u) {
      a <- exp(0.5 * (lgamma(1 / nu) - lgamma(3 / nu))); x <- (abs(u) / a)^nu
      g <- sign(u) * qnorm(log(0.5) + pgamma(x, 1 / nu, lower.tail = FALSE,
                                             log.p = TRUE),
                           lower.tail = FALSE, log.p = TRUE)
      c(g, exp(log(nu) - log(2 * a) - lgamma(1 / nu) - x - dnorm(g, log = TRUE)), 0) })
  ys <- log(y^2 + del) - pr$mu; Tn <- length(y)
  P <- sv^2 / (1 - phi^2); xi <- 0; u <- 0
  w <- numeric(Tn); Pp <- numeric(Tn)
  for (t in 1:Tn) {
    if (t == 1) {
      xp <- phi * xi; P_p <- phi^2 * P + sv^2
    } else {
      k <- gdv(u); A <- phi - 0.5 * sv * d * k[2] * u
      xp <- phi * xi + sv * d * k[1]
      P_p <- A^2 * P + sv^2 * ((1 - d^2) + d^2 * k[3])
    }
    Pp[t] <- P_p; f <- P_p + pr$sig_eps2; K <- P_p / f
    xi <- xp + K * (ys[t] - xp); P <- P_p - K * P_p; w[t] <- xi
    u <- (y[t] / pr$sigma_y) * exp(-0.5 * w[t])
  }
  list(w = w, Pp = Pp)
}

.lev_fit <- function(seed, errorType = "Gaussian", nu = NULL, rho = -0.6,
                     n = 800, sigv = 0.3) {
  set.seed(seed)
  y <- sim_svp(n, phi = 0.95, sigy = 1, sigv = sigv, errorType = errorType,
               nu = nu, leverage = TRUE, rho = rho)$y
  list(y = y, fit = suppressWarnings(svp(y, p = 1, errorType = errorType,
                                         leverage = TRUE)))
}

test_that("the proxy argument is gone from the four filter-consuming functions", {
  d <- .lev_fit(1, "Student-t", nu = 5)
  expect_silent(filter_svp(d$fit, method = "mixture"))
  expect_error(filter_svp(d$fit, proxy = "u"), "unused argument")
  expect_error(forecast_svp(d$fit, H = 2, proxy = "bayes_optimal"), "unused argument")
  expect_error(svp_IC(d$fit, proxy = "bayes_optimal"), "unused argument")
  expect_error(svp_AR_order(d$y, pmax = 2, proxy = "u"), "unused argument")
})

test_that("delta_p = 0 recovers the linear filter exactly (CKF, GMKF, smoother)", {
  for (et in c("Gaussian", "Student-t")) {
    d <- .lev_fit(11, et, nu = if (et == "Student-t") 5 else NULL)
    fit0 <- d$fit; fit0$rho <- 0                      # same parameters, no shift
    fitn <- suppressWarnings(svp(d$y, p = 1, errorType = et, leverage = FALSE))
    for (m in c("corrected", "mixture")) {
      a <- filter_svp(fit0, method = m); b <- filter_svp(fitn, method = m)
      for (f in c("w_filtered", "w_smoothed", "P_predicted", "P_filt_T", "loglik"))
        expect_identical(a[[f]], b[[f]], info = paste(et, m, f))
    }
  }
})

test_that("CKF leverage prediction matches the R reference of SVHT (48a)-(48b)", {
  cases <- list(list("Gaussian", NULL), list("Student-t", 5), list("GED", 1.5))
  for (cs in cases) {
    d <- .lev_fit(21, cs[[1]], nu = cs[[2]])
    ref <- .ckf_ref(d$fit, d$y)
    f <- filter_svp(d$fit, method = "corrected")
    expect_equal(f$w_filtered, ref$w, tolerance = 1e-10, info = cs[[1]])
    expect_equal(f$P_predicted, ref$Pp, tolerance = 1e-10, info = cs[[1]])
    expect_true(abs(d$fit$rho) > 0.1)                 # the shift is active
  }
})

test_that("h = 1 forecast equals the filter's own one-step prediction (CKF, GMKF)", {
  cases <- list(list("Gaussian", NULL), list("Student-t", 5), list("GED", 1.5))
  for (cs in cases) {
    d <- .lev_fit(31, cs[[1]], nu = cs[[2]])
    Tn <- length(d$y); fit1 <- d$fit; fit1$y <- d$y[-Tn]
    for (m in c("corrected", "mixture")) {
      fc <- forecast_svp(fit1, H = 1, filter_method = m)
      fl <- filter_svp(d$fit, method = m)
      expect_equal(fc$log_var_forecast[1], fl$w_predicted[Tn], tolerance = 1e-12,
                   info = paste(cs[[1]], m))
      expect_equal(fc$P_forecast[1], fl$P_predicted[Tn], tolerance = 1e-12,
                   info = paste(cs[[1]], m))
    }
  }
})

test_that("forecast covariance: (48a) at h = 1, marginal sigma_v^2 at h >= 2", {
  skip_on_cran()
  for (et in c("Gaussian", "Student-t")) {
    d <- .lev_fit(2026, et, nu = if (et == "Student-t") 5 else NULL, n = 1500, sigv = 0.4)
    fit <- d$fit; fc <- forecast_svp(fit, H = 5, filter_method = "corrected")
    phi <- as.numeric(fit$phi); sv <- fit$sigv; delta <- fit$rho
    uT <- fc$filter_output$zt[length(d$y)]
    if (et == "Gaussian") { dg <- 1; V <- 0 } else {
      nu <- fit$v; cc <- exp(lgamma((nu + 2) / 2) - lgamma((nu + 1) / 2)); s <- nu + uT^2
      dg <- sqrt(2) * nu * cc / s^1.5; V <- uT^2 * ((nu + 1) / s - 2 * cc^2 / s)
    }
    A1 <- phi - 0.5 * sv * delta * dg * uT
    q1 <- sv^2 * ((1 - delta^2) + delta^2 * V)
    P <- fc$P_forecast; P_T <- fc$filter_output$P_filt_T[1, 1]
    expect_equal(P[1] - A1^2 * P_T, q1, tolerance = 1e-8, info = et)
    expect_equal(P[2] - phi^2 * P[1], sv^2, tolerance = 1e-8, info = et)
    expect_equal(P[3] - phi^2 * P[2], sv^2, tolerance = 1e-8, info = et)
  }
})

test_that("leverage filters stay finite where 0.2.0.9005 diverged (rho at the bound)", {
  skip_on_cran()
  # Seeds on which the pre-0.2.0.9006 covariance (sigma_v^2 (1 - delta^2) alone)
  # sent the CKF and the GMKF to Inf: Gaussian SV(1), phi = 0.95, sigma_v = 0.5,
  # rho = -0.9, T = 2000, rho_hat at the truncation bound (-0.999 then, -0.99
  # since 0.2.0.9007).
  for (sd in c(6L, 19L, 23L)) {
    set.seed(sd)
    y <- sim_svp(2000, phi = 0.95, sigy = 1, sigv = 0.5, leverage = TRUE, rho = -0.9)$y
    fit <- suppressWarnings(svp(y, p = 1, leverage = TRUE))
    for (m in c("corrected", "mixture")) {
      f <- filter_svp(fit, method = m)
      expect_true(all(is.finite(f$w_filtered)) && max(abs(f$w_filtered)) < 20,
                  info = paste("seed", sd, m))
      expect_true(all(is.finite(f$w_smoothed)), info = paste("seed", sd, m))
      fc <- forecast_svp(fit, H = 22, filter_method = m)
      expect_true(all(is.finite(fc$log_var_forecast)), info = paste("seed", sd, m))
    }
  }
})

test_that("BPF Student-t recovery uses the posterior (Bug 3 fix)", {
  skip_on_cran()
  set.seed(3)
  y <- sim_svp(600, phi = 0.95, sigy = 1, sigv = 0.3,
                errorType = "Student-t", nu = 5,
                leverage = TRUE, rho = -0.4)$y
  fit <- suppressWarnings(svp(y, p = 1, errorType = "Student-t",
                               leverage = TRUE))
  # Just ensure BPF runs and returns finite log-likelihood after the
  # posterior-lambda fix.  A stronger empirical-mean test would require
  # exposing the recovery step internals.
  pf <- filter_svp(fit, method = "particle", M = 500, seed = 7)
  expect_true(is.finite(pf$loglik))
  expect_true(all(is.finite(pf$w_filtered)))
})

test_that("CKF returns a finite log-likelihood for every family, with and without leverage", {
  skip_on_cran()
  combos <- list(
    list(et = "Gaussian",  lev = FALSE),
    list(et = "Gaussian",  lev = TRUE,  rho = -0.3),
    list(et = "Student-t", lev = FALSE, nu = 5),
    list(et = "GED",       lev = FALSE, nu = 1.5)
  )
  for (i in seq_along(combos)) {
    set.seed(50 + i)
    args <- list(n = 500, phi = 0.95, sigy = 1, sigv = 0.3,
                 errorType = combos[[i]]$et, leverage = combos[[i]]$lev)
    if (!is.null(combos[[i]]$nu))  args$nu  <- combos[[i]]$nu
    if (!is.null(combos[[i]]$rho)) args$rho <- combos[[i]]$rho
    y <- do.call(sim_svp, args)$y
    fit <- suppressWarnings(svp(y, p = 1,
                                 errorType = combos[[i]]$et,
                                 leverage = combos[[i]]$lev))
    f <- filter_svp(fit, method = "corrected")
    expect_true(is.finite(f$loglik),
                info = paste("combo", i, combos[[i]]$et,
                              "leverage", combos[[i]]$lev))
  }
})
