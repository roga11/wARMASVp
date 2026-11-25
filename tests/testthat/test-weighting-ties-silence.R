# Weighting-matrix inverse, Monte Carlo p-value with ties, optimizer output.

test_that(".pinv inverts badly scaled matrices and is scale invariant", {
  pinv <- wARMASVp:::.pinv
  expect_equal(pinv(diag(c(1e9, 1))), diag(c(1e-9, 1)), tolerance = 1e-12)
  expect_equal(pinv(diag(c(1e-12, 1))), diag(c(1e12, 1)), tolerance = 1e-12)
  set.seed(1)
  X <- matrix(rnorm(60), 12, 5)
  W <- crossprod(X)
  expect_equal(pinv(W), solve(W), tolerance = 1e-10)
  D <- 10^c(-6, -2, 0, 3, 7)
  expect_equal(pinv(D * W * rep(D, each = 5)),
               (1 / D) * solve(W) * rep(1 / D, each = 5), tolerance = 1e-9)
})

test_that(".pinv warns on a singular matrix and rejects non-square input", {
  pinv <- wARMASVp:::.pinv
  S <- tcrossprod(c(1, 2, 3))
  expect_warning(G <- pinv(S), "numerically singular")
  expect_equal(S %*% G %*% S, S, tolerance = 1e-10)
  expect_error(pinv(matrix(1, 2, 3)), "square")
  expect_error(pinv(1:3), "square")
})

test_that("weighted GED leverage statistic does not depend on the scale of y", {
  set.seed(11)
  y <- sim_svp(1500, phi = 0.95, sigy = 1, sigv = 0.5, errorType = "GED",
               nu = 1.5, leverage = TRUE, rho = -0.5)$y
  s0 <- vapply(c(1e-4, 1, 1e3), function(cc) {
    set.seed(3)
    suppressWarnings(lmc_lev(cc * y, p = 1, J = 10, N = 9, errorType = "GED",
                             Amat = "Weighted")$s0)
  }, numeric(1))
  expect_true(all(is.finite(s0)))
  expect_equal(s0[1], s0[2], tolerance = 1e-8)
  expect_equal(s0[3], s0[2], tolerance = 1e-8)
})

test_that("simulated statistics tied with the observed one count toward the p-value", {
  # Gaussian data with nu_null at the estimator's upper bound: the observed
  # statistic is zero (reported at the 1e-10 floor) and so are the simulated
  # ones whose estimate is at the bound.
  set.seed(5)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.5)$y
  set.seed(6)
  out <- suppressWarnings(lmc_t(y, p = 1, J = 10, N = 19, nu_null = 500))
  expect_equal(out$s0, 1e-10)
  expect_true(any(out$sN == out$s0))
  expect_equal(out$pval, 1)
  expect_equal(out$pval, (1 + sum(out$sN >= out$s0)) / (19 + 1))
})

test_that("p-value without ties is (1 + #{sN >= s0}) / (N + 1)", {
  set.seed(7)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.5, errorType = "Student-t",
               nu = 5)$y
  set.seed(8)
  out <- suppressWarnings(lmc_t(y, p = 1, J = 10, N = 19, nu_null = 5))
  expect_equal(out$pval, (1 + sum(out$sN >= out$s0)) / 20)
  expect_equal(out$pval, (1 + sum(out$sN > out$s0)) / 20 +
                 sum(out$sN == out$s0) / 20)
})

# pso reports through message(), GenSA through cat(): count both.
n_printed <- function(expr) {
  n_msg <- 0L
  out <- capture.output(withCallingHandlers(
    invisible(suppressWarnings(expr)),
    message = function(m) {
      n_msg <<- n_msg + 1L
      invokeRestart("muffleMessage")
    }))
  length(out) + n_msg
}

test_that("silence = TRUE suppresses optimizer output without changing the result", {
  skip_if_not_installed("pso")
  set.seed(9)
  y <- sim_svp(600, phi = 0.9, sigy = 1, sigv = 0.5)$y
  run <- function(silence, method = "pso") {
    set.seed(10)
    res <- NULL
    n <- n_printed(res <- mmc_t(y, p = 1, J = 10, N = 9, nu_null = 8,
                                method = method, maxit = 6,
                                silence = silence))
    list(res = res, n = n)
  }
  quiet <- run(TRUE)
  loud <- run(FALSE)
  expect_identical(quiet$n, 0L)
  expect_gt(loud$n, 0L)
  expect_identical(quiet$res$pval, loud$res$pval)
  expect_identical(quiet$res$par, loud$res$par)
  expect_null(quiet$res$stats)
  expect_false(is.null(loud$res$stats))
  expect_error(mmc_t(y, nu_null = 8, N = 9, maxit = 2, silence = NA), "silence")
  skip_if_not_installed("GenSA")
  quiet_g <- run(TRUE, "GenSA")
  loud_g <- run(FALSE, "GenSA")
  expect_identical(quiet_g$n, 0L)
  expect_gt(loud_g$n, 0L)
  expect_identical(quiet_g$res$pval, loud_g$res$pval)
  expect_identical(quiet_g$res$par, loud_g$res$par)
})

test_that("silence is accepted by all four mmc_* functions", {
  skip_if_not_installed("pso")
  set.seed(12)
  yl <- sim_svp(600, phi = 0.9, sigy = 1, sigv = 0.5, leverage = TRUE,
                rho = -0.4)$y
  set.seed(13)
  expect_identical(n_printed(mmc_ar(yl, p_null = 1, p_alt = 2, N = 9,
                                    maxit = 4, silence = TRUE)), 0L)
  expect_identical(n_printed(mmc_lev(yl, N = 9, maxit = 4, silence = TRUE)), 0L)
  expect_identical(n_printed(mmc_ged(yl, N = 9, nu_null = 2, maxit = 4,
                                     silence = TRUE)), 0L)
  expect_gt(n_printed(mmc_lev(yl, N = 9, maxit = 4)), 0L)
})

test_that("mmc_lev default search region: fixed nu half-width, sigma_y rule", {
  skip_if_not_installed("pso")
  run <- function(y, eT, eps, Amat = NULL) {
    set.seed(21)
    suppressWarnings(mmc_lev(y, N = 9, errorType = eT, eps = eps, Amat = Amat,
                             maxit = 8, silence = TRUE))
  }
  same <- function(a, b) {
    expect_identical(a$pval, b$pval)
    expect_identical(a$par, b$par)
  }
  set.seed(20)
  yt <- sim_svp(800, phi = 0.9, sigy = 1, sigv = 0.5, errorType = "Student-t",
                nu = 20, leverage = TRUE, rho = -0.4)$y
  ft <- suppressWarnings(svp(yt, 1, 10, leverage = TRUE, errorType = "Student-t"))
  d <- run(yt, "Student-t", NULL)
  same(d, run(yt, "Student-t", c(0.3, 0, 0.3, 0.3)))
  same(run(yt, "Student-t", c(0.2, 0, 0.2)),
       run(yt, "Student-t", c(0.2, 0, 0.2, 0.3)))
  expect_lte(abs(d$par[4] - ft$v), 0.3 + 1e-12)
  expect_identical(d$par[2], ft$sigy)

  set.seed(22)
  yg <- sim_svp(800, phi = 0.9, sigy = 1, sigv = 0.5, errorType = "GED",
                nu = 1.5, leverage = TRUE, rho = -0.4)$y
  fg <- suppressWarnings(svp(yg, 1, 10, leverage = TRUE, errorType = "GED"))
  # identity weighting: sigma_y is searched
  same(run(yg, "GED", NULL), run(yg, "GED", c(0.3, 0.3, 0.3, 0.3)))
  # HAC weighting: sigma_y is held at its estimate
  h <- run(yg, "GED", NULL, "Weighted")
  same(h, run(yg, "GED", c(0.3, 0, 0.3, 0.3), "Weighted"))
  expect_identical(h$par[2], fg$sigy)
  # a user-supplied zero is respected
  expect_identical(run(yg, "GED", c(0.3, 0, 0.3, 0.3))$par[2], fg$sigy)
})

test_that("Kendall leverage test refuses HAC weighting, allows the others", {
  set.seed(11)
  y <- sim_svp(400, phi = 0.9, sigy = 1, sigv = 1, leverage = TRUE, rho = -0.5)$y
  msg <- "not supported with HAC weighting"
  expect_error(lmc_lev(y, N = 9, rho_type = "kendall", Amat = "Weighted"), msg)
  expect_error(lmc_lev(y, N = 9, rho_type = "kendall", Bartlett = TRUE), msg)
  expect_error(mmc_lev(y, N = 9, rho_type = "kendall", Amat = "Weighted"), msg)
  expect_error(mmc_lev(y, N = 9, rho_type = "kendall", Bartlett = TRUE), msg)
  expect_s3_class(suppressWarnings(lmc_lev(y, N = 9, rho_type = "kendall")), "svp_test")
  expect_s3_class(suppressWarnings(
    lmc_lev(y, N = 9, rho_type = "kendall", Amat = diag(4))), "svp_test")
  expect_s3_class(suppressWarnings(
    lmc_lev(y, N = 9, rho_type = "pearson", Amat = "Weighted")), "svp_test")
})
