# del = 0 default (log(y^2) with no offset), the zero-return guard, and the
# fit-carried del used by filter_svp() / forecast_svp() / svp_IC().

.del_y <- function(seed = 11, n = 600) {
  set.seed(seed)
  sim_svp(n, phi = 0.95, sigy = 1, sigv = 0.5)$y
}

test_that("del = 0 is the default and equals the del -> 0 limit", {
  y <- .del_y()
  f0 <- svp(y, p = 1)
  expect_identical(f0$del, 0)
  expect_identical(f0$theta, svp(y, p = 1, del = 1e-300)$theta)
  for (et in c("Student-t", "GED")) {
    a <- suppressWarnings(svp(y, p = 1, errorType = et))
    b <- suppressWarnings(svp(y, p = 1, errorType = et, del = 1e-300))
    expect_identical(a$theta, b$theta, info = et)
  }
})

test_that("default estimates are invariant to the scale of y", {
  y <- .del_y()
  for (et in c("Gaussian", "Student-t", "GED")) {
    a <- suppressWarnings(svp(y, p = 1, errorType = et))
    b <- suppressWarnings(svp(1e-4 * y, p = 1, errorType = et))
    expect_equal(b$phi, a$phi, tolerance = 1e-8, info = et)
    expect_equal(b$sigv, a$sigv, tolerance = 1e-8, info = et)
    expect_equal(b$sigy, 1e-4 * a$sigy, tolerance = 1e-8, info = et)
  }
})

test_that("zeros in y with del = 0 give an informative error everywhere", {
  y0 <- .del_y(); y0[c(10, 300)] <- 0
  msg <- "contains 2 zero value\\(s\\) \\(first at index 10\\).*Set 'del'"
  expect_error(svp(y0), msg)
  expect_error(svp(y0, errorType = "Student-t"), msg)
  expect_error(svp(y0, errorType = "GED", leverage = TRUE), msg)
  expect_error(svp_AR_order(y0, pmax = 2), msg)
  expect_error(lmc_ar(y0, p_null = 1, p_alt = 2, N = 9), msg)
  expect_error(lmc_lev(y0, N = 9), msg)
  expect_error(lmc_t(y0, nu_null = 5, N = 9), msg)
  expect_error(lmc_ged(y0, nu_null = 1.5, N = 9), msg)
  expect_error(mmc_ar(y0, p_null = 1, p_alt = 2, N = 9), msg)
  expect_error(mmc_lev(y0, N = 9), msg)
  expect_error(mmc_t(y0, nu_null = 5, N = 9), msg)
  expect_error(mmc_ged(y0, nu_null = 1.5, N = 9), msg)
})

test_that("a positive del handles zeros and is carried by the fit", {
  y0 <- .del_y(); y0[c(10, 300)] <- 0
  fit <- svp(y0, p = 1, del = 1e-8)
  expect_identical(fit$del, 1e-8)
  f <- filter_svp(fit)
  expect_true(all(is.finite(f$w_filtered)))
  expect_identical(f$loglik, filter_svp(fit, del = 1e-8)$loglik)
  expect_true(all(is.finite(forecast_svp(fit, H = 2)$log_var_forecast)))
  expect_true(all(is.finite(svp_IC(fit))))
  expect_true(all(is.finite(filter_svp(fit, method = "mixture")$w_filtered)))
  expect_true(all(is.finite(filter_svp(fit, method = "particle", M = 200)$w_filtered)))
  # an explicit del = 0 on such a fit is refused
  expect_error(filter_svp(fit, del = 0), "zero value")
  expect_error(forecast_svp(fit, del = 0), "zero value")
  expect_error(svp_IC(fit, del = 0), "zero value")
  ord <- svp_AR_order(y0, pmax = 2, del = 1e-8)
  expect_true(all(is.finite(ord$IC)))
})

test_that("del must be a non-negative number", {
  y <- .del_y()
  expect_error(svp(y, del = -1e-10), "non-negative")
  expect_error(svp(y, del = NA_real_), "non-negative")
  expect_error(svp(y, del = c(0, 1)), "non-negative")
  expect_error(lmc_lev(y, N = 9, del = -1), "non-negative")
  expect_error(filter_svp(svp(y), del = -1), "non-negative")
})
