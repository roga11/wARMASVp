test_that("svp_IC returns the 4 default criteria for Gaussian SV(1) without leverage", {
  set.seed(101)
  y <- sim_svp(800, phi = 0.95, sigy = 1, sigv = 0.4)$y
  fit <- svp(y, p = 1)
  ic <- svp_IC(fit)
  expect_named(ic, c("BIC_Kalman", "AIC_Kalman", "BIC_HR", "AIC_HR"),
               ignore.order = TRUE)
  expect_true(all(is.finite(ic)))
})

test_that("svp_IC handles all errorType x leverage combinations", {
  combos <- expand.grid(
    errorType = c("Gaussian", "Student-t", "GED"),
    leverage = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(combos))) {
    set.seed(200 + i)
    args <- list(n = 800, phi = 0.95, sigy = 1, sigv = 0.4,
                 errorType = combos$errorType[i],
                 leverage = combos$leverage[i])
    if (combos$errorType[i] == "Student-t") args$nu <- 5
    if (combos$errorType[i] == "GED")       args$nu <- 1.5
    if (combos$leverage[i])                 args$rho <- -0.4
    y <- do.call(sim_svp, args)$y
    fit <- suppressWarnings(svp(y, p = 1,
                                 errorType = combos$errorType[i],
                                 leverage = combos$leverage[i]))
    ic <- svp_IC(fit)
    expect_length(ic, 4L)
    # at least one finite criterion
    expect_true(any(is.finite(ic)),
                info = sprintf("Row %d: %s, leverage=%s",
                               i, combos$errorType[i], combos$leverage[i]))
  }
})

test_that("svp_IC criteria subset selection works", {
  set.seed(102)
  y <- sim_svp(800, phi = 0.95, sigy = 1, sigv = 0.4)$y
  fit <- svp(y, p = 1)

  ic_one <- svp_IC(fit, criteria = "BIC_Kalman")
  expect_named(ic_one, "BIC_Kalman")
  expect_length(ic_one, 1L)

  ic_two <- svp_IC(fit, criteria = c("BIC_YW", "AIC_YW"))
  expect_named(ic_two, c("BIC_YW", "AIC_YW"), ignore.order = TRUE)
  expect_length(ic_two, 2L)
})

test_that("svp_IC computes all 8 criteria when explicitly requested", {
  set.seed(106)
  y <- sim_svp(800, phi = 0.95, sigy = 1, sigv = 0.4)$y
  fit <- svp(y, p = 1)
  all_crit <- c("BIC_Kalman", "AIC_Kalman", "AICc_Kalman",
                "BIC_Whittle",
                "BIC_HR", "AIC_HR",
                "BIC_YW", "AIC_YW")
  ic_all <- svp_IC(fit, criteria = all_crit)
  expect_length(ic_all, 8L)
  expect_named(ic_all, all_crit, ignore.order = TRUE)
  # at least one finite criterion in each opt-in family
  expect_true(is.finite(ic_all["AICc_Kalman"]))
  expect_true(is.finite(ic_all["BIC_Whittle"]))
  expect_true(is.finite(ic_all["BIC_YW"]))
  expect_true(is.finite(ic_all["AIC_YW"]))
})

test_that("svp_IC errors on non-svp object", {
  expect_error(svp_IC(list(phi = 0.5)),
               "must be an svp")
})

test_that("svp_AR_order at SV(1) DGP selects p=1 by Kalman criteria", {
  skip_on_cran()
  set.seed(103)
  y <- sim_svp(2000, phi = 0.95, sigy = 1, sigv = 0.5)$y
  res <- svp_AR_order(y, pmax = 3L)
  expect_true(is.matrix(res$IC))
  expect_equal(dim(res$IC), c(4L, 3L))
  expect_named(res$argmin,
               c("BIC_Kalman", "AIC_Kalman", "BIC_HR", "AIC_HR"),
               ignore.order = TRUE)
  # Kalman criteria should select p=1 for SV(1) truth
  expect_equal(unname(res$argmin["BIC_Kalman"]), 1L)
  expect_equal(unname(res$argmin["AIC_Kalman"]), 1L)
})

test_that("svp_AR_order pmax=1 returns 1-column IC matrix", {
  set.seed(104)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.4)$y
  res <- svp_AR_order(y, pmax = 1L)
  expect_equal(ncol(res$IC), 1L)
  expect_equal(length(res$fits), 1L)
  expect_true(all(res$argmin == 1L | is.na(res$argmin)))
})

test_that("internal helpers produce sensible values", {
  set.seed(105)
  y <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.5)$y
  fit <- svp(y, p = 1)

  # n_params: SV(1) Gaussian no-lev = 1 + 2 = 3
  expect_equal(wARMASVp:::.svp_n_params(1L, "Gaussian", FALSE), 3L)
  # SV(2) Student-t with leverage = 2 + 2 + 1 + 1 = 6
  expect_equal(wARMASVp:::.svp_n_params(2L, "Student-t", TRUE), 6L)

  # sigma_eps^2 for Gaussian
  expect_equal(wARMASVp:::.svp_sigma_eps2("Gaussian"), pi^2 / 2)
  # for Student-t at nu=5
  expect_equal(wARMASVp:::.svp_sigma_eps2("Student-t", nu = 5),
               trigamma(0.5) + trigamma(2.5))

  # acov returns max_lag+1 finite values
  acov <- wARMASVp:::.svp_acov_ystar(y, max_lag = 3L)
  expect_length(acov, 4L)
  expect_true(all(is.finite(acov)))
  expect_true(acov[1] > 0)  # variance positive

  # residual variance is finite and positive for a real fit
  rv <- wARMASVp:::.svp_residual_var(fit, y)
  expect_true(is.finite(rv) && rv > 0)
})


test_that("Hannan-Rissanen helpers run and produce sensible variance", {
  skip_on_cran()
  set.seed(1234)
  y <- sim_svp(1500, phi = 0.95, sigy = 1, sigv = 0.5)$y
  fit <- svp(y, p = 1)

  # Direct call to the wrapper
  hr <- wARMASVp:::.svp_hr_residual_var(fit, y)
  expect_true(is.finite(hr$sigma2))
  expect_true(hr$sigma2 > 0)
  expect_true(hr$T_eff > 100L)

  # Stage 1 alone: should produce residuals
  ystar <- log(y^2 + 1e-10); ystar <- ystar - mean(ystar)
  s1 <- wARMASVp:::.svp_hr_stage1(ystar)
  expect_true(!is.null(s1$eps))
  expect_true(s1$L > 0L)
  expect_true(length(s1$eps) > 0L)
})

test_that("BIC_HR / AIC_HR appear in svp_IC output and are finite", {
  skip_on_cran()
  set.seed(2345)
  y <- sim_svp(1500, phi = 0.95, sigy = 1, sigv = 0.5)$y
  fit <- svp(y, p = 2)
  ic <- svp_IC(fit, criteria = c("BIC_HR", "AIC_HR"))
  expect_named(ic, c("BIC_HR", "AIC_HR"), ignore.order = TRUE)
  expect_true(all(is.finite(ic)))
})
