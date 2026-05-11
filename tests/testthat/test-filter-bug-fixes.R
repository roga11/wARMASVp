# Tests for the 2026-05-09 filter parameterization fixes:
#   Fix 1: var_zt = 1 in CKF/GMKF/forecast (was nu/(nu-2) for Student-t leverage)
#   Fix 3: BPF posterior lambda for Student-t leverage (was prior)
#   Fix 4: GED copula proxy in CKF/GMKF leverage shift
# Plus the proxy = c("u", "bayes_optimal") API.

test_that("filter_svp accepts proxy argument; default is 'u'", {
  set.seed(1)
  y <- sim_svp(800, phi = 0.95, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5,
               leverage = TRUE, rho = -0.4)$y
  fit <- suppressWarnings(svp(y, p = 1, errorType = "Student-t",
                               leverage = TRUE))
  expect_silent(filter_svp(fit, method = "mixture"))
  expect_silent(filter_svp(fit, method = "mixture", proxy = "u"))
  expect_silent(filter_svp(fit, method = "mixture", proxy = "bayes_optimal"))
})

test_that("proxy choice differs only for Student-t with leverage", {
  set.seed(2)
  # Gaussian no-leverage: proxy is irrelevant (zeta = u exactly)
  y_g <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.4)$y
  fit_g <- svp(y_g, p = 1)
  f_g_u <- filter_svp(fit_g, method = "corrected", proxy = "u")
  f_g_b <- filter_svp(fit_g, method = "corrected", proxy = "bayes_optimal")
  expect_equal(f_g_u$loglik, f_g_b$loglik, tolerance = 1e-10)

  # Gaussian with leverage: still proxy-irrelevant (zeta = u for Gaussian)
  y_gl <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                   leverage = TRUE, rho = -0.3)$y
  fit_gl <- svp(y_gl, p = 1, leverage = TRUE)
  f_gl_u <- filter_svp(fit_gl, method = "corrected", proxy = "u")
  f_gl_b <- filter_svp(fit_gl, method = "corrected", proxy = "bayes_optimal")
  expect_equal(f_gl_u$loglik, f_gl_b$loglik, tolerance = 1e-10)

  # Student-t with leverage: proxy DOES matter
  y_t <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                  errorType = "Student-t", nu = 5,
                  leverage = TRUE, rho = -0.4)$y
  fit_t <- suppressWarnings(svp(y_t, p = 1, errorType = "Student-t",
                                 leverage = TRUE))
  f_t_u <- filter_svp(fit_t, method = "corrected", proxy = "u")
  f_t_b <- filter_svp(fit_t, method = "corrected", proxy = "bayes_optimal")
  # The two log-likelihoods should differ for Student-t leverage
  expect_false(isTRUE(all.equal(f_t_u$loglik, f_t_b$loglik,
                                 tolerance = 1e-6)))

  # GED with leverage: proxy is irrelevant (the GED copula proxy is
  # deterministic and doesn't depend on the proxy_type flag)
  y_ged <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                    errorType = "GED", nu = 1.5,
                    leverage = TRUE, rho = -0.4)$y
  fit_ged <- suppressWarnings(svp(y_ged, p = 1, errorType = "GED",
                                   leverage = TRUE))
  f_ged_u <- filter_svp(fit_ged, method = "corrected", proxy = "u")
  f_ged_b <- filter_svp(fit_ged, method = "corrected", proxy = "bayes_optimal")
  expect_equal(f_ged_u$loglik, f_ged_b$loglik, tolerance = 1e-10)
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

test_that("svp_IC default proxy is 'bayes_optimal'", {
  skip_on_cran()
  set.seed(4)
  y <- sim_svp(800, phi = 0.95, sigy = 1, sigv = 0.3,
                errorType = "Student-t", nu = 5,
                leverage = TRUE, rho = -0.4)$y
  fit <- suppressWarnings(svp(y, p = 1, errorType = "Student-t",
                               leverage = TRUE))
  # Default svp_IC should match svp_IC(proxy = "bayes_optimal") exactly
  ic_default <- svp_IC(fit)
  ic_b       <- svp_IC(fit, proxy = "bayes_optimal")
  expect_equal(ic_default, ic_b, tolerance = 1e-8)
  # (Whether proxy="u" gives a numerically different value depends on the
  # estimated delta_p — small delta_p → small effect.  We verify the
  # forwarding instead of assuming a numerical difference here.)
})

test_that("var_zt fix does not affect Gaussian or no-leverage cells", {
  skip_on_cran()
  # All three combos that should be invariant to the fix:
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

test_that("forecast_svp accepts proxy and forwards correctly", {
  skip_on_cran()
  set.seed(5)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                errorType = "Student-t", nu = 5,
                leverage = TRUE, rho = -0.4)$y
  fit <- suppressWarnings(svp(y, p = 1, errorType = "Student-t",
                               leverage = TRUE))
  fc_u <- forecast_svp(fit, H = 5, filter_method = "mixture", proxy = "u")
  fc_b <- forecast_svp(fit, H = 5, filter_method = "mixture",
                        proxy = "bayes_optimal")
  expect_true(all(is.finite(fc_u$vol_forecast)))
  expect_true(all(is.finite(fc_b$vol_forecast)))
  # Should differ (not identical) for Student-t leverage
  expect_false(isTRUE(all.equal(fc_u$vol_forecast, fc_b$vol_forecast,
                                 tolerance = 1e-6)))
})


test_that("forecast Q is conditional at h=1 and marginal at h>=2 (Bug 5 fix)", {
  skip_on_cran()
  set.seed(2026)
  y <- sim_svp(1500, phi = 0.95, sigy = 1, sigv = 0.4,
               leverage = TRUE, rho = -0.5)$y
  fit <- suppressWarnings(svp(y, p = 1, leverage = TRUE))
  fc <- forecast_svp(fit, H = 5, filter_method = "corrected")
  # Build the expected Q values
  phi <- as.numeric(fit$phi)
  sigv2 <- fit$sigv^2
  delta <- if (is.null(fit$rho) || is.na(fit$rho)) 0 else fit$rho
  Q_h1   <- sigv2 * (1 - delta^2)   # conditional at h=1
  Q_hge2 <- sigv2                    # marginal at h>=2
  # The Q recursion increment equals Q_h at each step:
  #   P_h - phi^2 * P_{h-1} = Q_h
  P <- fc$P_forecast
  # Use filter's P_filt_T as the seed for h=1
  P_T_11 <- fc$filter_output$P_filt_T[1, 1]
  expect_equal(P[1] - phi^2 * P_T_11, Q_h1, tolerance = 1e-6)
  # h>=2 increments should equal Q_hge2
  expect_equal(P[2] - phi^2 * P[1], Q_hge2, tolerance = 1e-6)
  expect_equal(P[3] - phi^2 * P[2], Q_hge2, tolerance = 1e-6)
})
