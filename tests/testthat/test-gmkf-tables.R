# GMKF mixture-table tests (2026-08-05): the precomputed "optimal" tables,
# their interpolation, the mixture-source dispatch in filter_svp(), and
# backwards compatibility of the legacy "em"/"ksc" paths.

# ---- 1. shipped tables are byte-identical to the built artifacts ---------- #

test_that("shipped table files match their build md5s", {
  expected <- c(gaussian = "b49f28c1afebb625298779df525134fd",
                ged      = "a3063cf6a4421a911f0c53691a57cd31",
                studentt = "cd3cbe09ee274ba689eb99b65dc34db4")
  for (fam in names(expected)) {
    path <- system.file("extdata", paste0(fam, "_mixtures.rds"),
                        package = "wARMASVp")
    expect_true(nzchar(path), info = fam)
    expect_identical(unname(tools::md5sum(path)), unname(expected[fam]),
                     info = fam)
  }
})

# ---- 2. interpolation reproduces stored nodes exactly --------------------- #

test_that("gmkf mixture at a stored node equals the stored fit", {
  for (fam in c("studentt", "ged")) {
    tab <- gmkf_table(fam)
    itp <- wARMASVp:::gmkf_interpolant(tab, fam)
    tb  <- tab[["7"]]
    i   <- length(tb$nu) %/% 2L      # an interior node
    g   <- wARMASVp:::gmkf_mixture(itp, tb$nu[i], 7)
    expect_equal(g$w, unname(tb$w[i, ]), tolerance = 1e-12, info = fam)
    expect_equal(g$m, unname(tb$m[i, ]), tolerance = 1e-12, info = fam)
    expect_equal(g$v, unname(tb$v[i, ]), tolerance = 1e-12, info = fam)
  }
})

# ---- 3. proper mixture at arbitrary nu, incl. the estimator clamps -------- #

test_that("interpolated mixtures are proper at arbitrary nu and at the clamps", {
  cases <- list(studentt = c(2.01, 2.81, 17.3, 500),
                ged      = c(0.1, 0.63, 1.22, 20))
  for (fam in names(cases)) {
    itp <- wARMASVp:::gmkf_interpolant(gmkf_table(fam), fam)
    for (nu in cases[[fam]]) {
      for (K in c(7, 10, 12, 16)) {
        g <- wARMASVp:::gmkf_mixture(itp, nu, K)
        expect_length(g$w, K)
        expect_true(all(g$w > 0) && abs(sum(g$w) - 1) < 1e-12,
                    info = sprintf("%s nu=%g K=%d", fam, nu, K))
        expect_true(all(g$v > 0) && all(is.finite(g$m)),
                    info = sprintf("%s nu=%g K=%d", fam, nu, K))
      }
    }
  }
})

# ---- fits shared by the dispatch tests ------------------------------------ #

.gmkf_fit <- function(seed, errorType = "Gaussian", nu = NULL, leverage = FALSE) {
  set.seed(seed)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.5, errorType = errorType,
               nu = nu, leverage = leverage, rho = if (leverage) -0.5 else 0)$y
  suppressWarnings(svp(y, p = 1, errorType = errorType, leverage = leverage))
}
fitA <- .gmkf_fit(101)                                     # Gaussian
fitB <- .gmkf_fit(302, "Student-t", nu = 3)                # Student-t, interior nu
fitC <- .gmkf_fit(103, "GED", nu = 1.5)                    # GED
fitD <- .gmkf_fit(104, "Student-t", nu = 3, leverage = TRUE)

# ---- 4. mixture = "em" reproduces 0.2.0.9002 exactly ---------------------- #

test_that("mixture = 'em' matches frozen 0.2.0.9002 logliks", {
  # Literals frozen from the installed 0.2.0.9002 on 2026-08-05 (scratch
  # script freeze_baselines_0202_9002.R), same seeds/fits as above.
  expect_equal(filter_svp(fitA, method = "mixture", mixture = "em",
                          K = 7, seed = 42)$loglik,
               -670.756399171503, tolerance = 1e-10)
  expect_equal(filter_svp(fitA, method = "mixture", mixture = "em",
                          K = 5, seed = 42)$loglik,
               -667.752180421772, tolerance = 1e-10)   # gaussian EM at K != 7
  expect_equal(filter_svp(fitB, method = "mixture", mixture = "em",
                          K = 7, seed = 42)$loglik,
               -722.224040723902, tolerance = 1e-10)
  expect_equal(filter_svp(fitB, method = "mixture", mixture = "em",
                          K = 3, seed = 42)$loglik,
               -721.941900006826, tolerance = 1e-10)   # arbitrary-K EM
  expect_equal(filter_svp(fitC, method = "mixture", mixture = "em",
                          K = 7, seed = 42)$loglik,
               -669.786206686606, tolerance = 1e-10)
  expect_equal(filter_svp(fitD, method = "mixture", mixture = "em",
                          K = 7, seed = 42)$loglik,
               -683.571674389611, tolerance = 1e-10)   # Student-t leverage
})

test_that("mixture = 'em' IC and forecast match frozen 0.2.0.9002 values", {
  expect_equal(unname(svp_IC(fitC, criteria = c("BIC_Kalman", "AIC_Kalman"),
                             mixture = "em")),
               c(1364.430845766901, 1347.572413373212), tolerance = 1e-10)
  expect_equal(forecast_svp(fitD, H = 2, filter_method = "mixture",
                            mixture = "em")$log_var_forecast,
               c(-1.644374704648, -1.445872235208), tolerance = 1e-10)
})

# ---- 5. "ksc" coincides with "em" for Gaussian K = 7; errors otherwise ---- #

test_that("mixture = 'ksc' is the Gaussian K = 7 table and nothing else", {
  f_ksc <- filter_svp(fitA, method = "mixture", mixture = "ksc", K = 7)
  f_em  <- filter_svp(fitA, method = "mixture", mixture = "em", K = 7)
  expect_identical(f_ksc$loglik, f_em$loglik)
  expect_identical(f_ksc$mixture, f_em$mixture)
  expect_error(filter_svp(fitB, method = "mixture", mixture = "ksc", K = 7),
               "ksc")
  expect_error(filter_svp(fitA, method = "mixture", mixture = "ksc", K = 10),
               "ksc")
})

# ---- 6. "optimal" is deterministic -------------------------------------- #

test_that("mixture = 'optimal' is invariant to seed", {
  l1 <- filter_svp(fitB, method = "mixture", seed = 1)$loglik
  l2 <- filter_svp(fitB, method = "mixture", seed = 2)$loglik
  expect_identical(l1, l2)
  expect_true(is.finite(l1))
})

# ---- 7. K validation: curated error, also through svp_IC ------------------ #

test_that("optimal path errors on unsupported K, incl. through svp_IC", {
  expect_error(filter_svp(fitA, method = "mixture", K = 9),
               "7, 10, 12, 16")
  expect_error(filter_svp(fitB, method = "mixture", K = 7.5),
               "7, 10, 12, 16")
  # must ERROR through svp_IC (not silently return NA Kalman ICs)
  expect_error(svp_IC(fitB, criteria = "BIC_Kalman", K = 9),
               "7, 10, 12, 16")
  # but non-mixture filter methods with odd K keep working
  ic <- svp_IC(fitB, criteria = "BIC_Kalman", filter_method = "corrected",
               K = 9)
  expect_true(is.finite(ic["BIC_Kalman"]))
})

# ---- 8. forwarding through forecast_svp / svp_IC / svp_AR_order ----------- #

test_that("mixture forwards through forecast_svp, svp_IC and svp_AR_order", {
  ic_opt <- svp_IC(fitB, criteria = "BIC_Kalman")
  ic_em  <- svp_IC(fitB, criteria = "BIC_Kalman", mixture = "em")
  expect_false(identical(ic_opt, ic_em))
  expect_identical(ic_opt, svp_IC(fitB, criteria = "BIC_Kalman", seed = 7))

  fc_opt <- forecast_svp(fitB, H = 2, filter_method = "mixture")
  fc_em  <- forecast_svp(fitB, H = 2, filter_method = "mixture",
                         mixture = "em")
  expect_false(identical(fc_opt$log_var_forecast, fc_em$log_var_forecast))

  set.seed(101)
  y <- sim_svp(400, phi = 0.9, sigy = 1, sigv = 0.5)$y
  r1 <- svp_AR_order(y, pmax = 2L, criteria = "BIC_Kalman", seed = 1)
  r2 <- svp_AR_order(y, pmax = 2L, criteria = "BIC_Kalman", seed = 2)
  expect_identical(r1$IC, r2$IC)   # optimal default: seed-invariant sweep
})

# ---- 9. RNG hygiene of the optimal path ----------------------------------- #

test_that("optimal path preserves a pre-existing .Random.seed exactly", {
  # Value-preservation given a pre-existing seed. (Seed-absence cannot be
  # asserted through filter_svp: Rcpp's RNGScope re-materialises .Random.seed
  # on any C++ kernel return -- see the note in test-rng-hygiene.R.)
  set.seed(2718)
  s_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(filter_svp(fitB, method = "mixture", K = 7)$loglik)
  s_after <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  expect_identical(s_before, s_after)
})

# ---- 10. accessor + output metadata ---------------------------------------- #

test_that("gmkf_table() returns provenance; filter output records the source", {
  for (fam in c("gaussian", "ged", "studentt")) {
    prov <- attr(gmkf_table(fam), "provenance")
    expect_true(is.list(prov) && all(c("script_md5", "helpers_md5") %in%
                                       names(prov)), info = fam)
  }
  f_opt <- filter_svp(fitB, method = "mixture")
  expect_identical(f_opt$mixture_source, "optimal")
  expect_identical(f_opt$mixture$distribution, "student_t")  # package code
  expect_null(f_opt$mixture$KL_div)
  f_em <- filter_svp(fitB, method = "mixture", mixture = "em")
  expect_identical(f_em$mixture_source, "em")
  expect_identical(f_em$mixture$distribution, "student_t")
  f_ckf <- filter_svp(fitB, method = "corrected")
  expect_null(f_ckf$mixture_source)
})
