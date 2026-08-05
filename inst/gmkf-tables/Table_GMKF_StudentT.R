## ===========================================================================
## Table_GMKF_StudentT.R
##
## Builds KL-optimal Gaussian-mixture approximations to the density of the
## centred log squared innovation
##
##     x = log(z_t^2) - E[log(z_t^2)],     z_t ~ Student-t(nu),
##
## for K = 7, 10, 12 and 16 components, tabulated over a grid of nu.  This is
## the measurement-error density of the linearised SV(p) observation equation
##
##     y*_t = log(y_t^2) - mu = h_t + x_t,   mu = log(sigma_y^2) + E[log(z_t^2)],
##
## and the mixture is what the Gaussian-mixture Kalman filter (GMKF) uses in
## place of it.
##
## CENTRING.  Everything is on the CENTRED scale, matching y*_t above and the
## mixtures the package fits for the other error families.  The innovation is
## the UNSTANDARDISED Student-t, so Var(z) = nu/(nu-2), and z^2 ~ F(1, nu) with
##     E[log(z_t^2)] = psi(1/2) - psi(nu/2) + log(nu).
##
## NU GRID.  nu in [2.01, 500], 500 nodes spaced uniformly in 1/nu.  The range
## is the estimator's own clamp (wARMASVp/src/rootfinding.cpp: nu_lower = 2.01,
## nu_upper = 500.0; the same bounds appear in R/svp_internal.R).  A table over
## a narrower range would leave reachable estimates untabulated.
##
## Why 1/nu.  What matters for continuation is that consecutive targets are
## about equally far apart, measured as a distance between densities rather
## than as a gap in nu.  On this family, spacing uniformly in 1/nu gives a
## worst-to-median Hellinger step of 1.4; uniform in log(nu) gives 9.2 and
## uniform in nu gives 6243.  It also concentrates nodes where the density
## actually changes: 401 of the 500 sit at nu <= 10, and only 9 above nu = 100,
## because beyond that the density is already close to its nu -> infinity
## limit and successive nodes barely differ.
## The grid is DEFINED in 1/nu and inverted; everything stored, printed and
## returned is plain nu.
##
## TWO GRIDS, and they are different things:
##   QUADRATURE POINTS  values of x at which the target density is evaluated to
##                      form the objective.  2000 per fit, rebuilt at each nu
##                      because the target changes with nu.
##   NODES              the tabulated values of nu.  200 of them here.
##
## METHOD
##   1. Deterministic quadrature.  At each nu the target is evaluated at 2000
##      quadrature points spanning all but ~1e-13 of its mass, giving weights
##      p_i, and the objective is
##          F(w,m,v) = - sum_i p_i log g(x_i),        F = KL + constant.
##      The mixture is fitted to a DENSITY, not to data: the quadrature points
##      and their weights play the role of a sample.  No sampling and no RNG
##      anywhere, so the result is bit-reproducible and carries no seed.  This
##      is a one-off table-building computation; nothing at filter runtime
##      repeats it.
##   2. Unconstrained gauge-fixed chart.  w = softmax(eta) with eta_1 = 0 held
##      fixed, v = exp(ell).  Gradient and Hessian are analytic.
##   3. Weighted EM, then a trust-region Newton polish.  EM is used only as a
##      globaliser: monotone and unable to diverge, so it carries a crude start
##      into a good basin, but its rate approaches 1 and it stalls well short
##      of the optimum.  The trust-region step converges it, solving each
##      subproblem by exact eigendecomposition so indefinite Hessians and the
##      hard case are handled rather than damped away.
##   4. Continuation in nu.  The mixture at the largest nu is built from
##      scratch by a nested ladder in K: K = 2, 3, ..., 16, where each rung
##      starts by splitting each component of the polished K-1 fit in turn and
##      keeps whichever candidate polishes to the lowest objective.  Every
##      other node is then reached by walking DOWN the grid, each node
##      warm-started from its neighbour.
##      The walk starts at the largest nu because that is where the mixture is
##      most tightly determined: log(z^2) is hardest to approximate there, so
##      every component carries real weight.  At the small-nu end the objective
##      has many more local minima, so a chain arriving there with a good fit in
##      hand is on safer ground than one starting there.
##   4b. Optional check at the largest nu.  As nu grows, log(z^2) tends to
##      log(chi^2_1), the density Table_GMKF_Gaussian.R solves by a separate
##      route.  The limit is not reached at nu = 500 -- the two densities still
##      differ by about 3.5e-04 at their widest -- so this is a tolerance check
##      rather than an identity: the mixtures should be close, and a large
##      disagreement would indicate a problem.  It is skipped if that file is
##      absent, and it never alters the table.
##   5. Component tracking.  After each nu step the new components are matched
##      to the previous node's by optimal assignment, so component k follows a
##      continuous path along the grid.  Sorting by mean instead would swap
##      identities wherever two component means cross.
##   6. Adaptive quadrature refinement.  A variance floor vmin = (2*Delta)^2 is
##      imposed, Delta being the quadrature spacing, because the DISCRETISED
##      objective is unbounded: a component whose variance shrinks to zero
##      while centred on a quadrature point drives F to -infinity, which is an
##      artefact of the discretisation.  If a fit lands ON the floor the grid
##      is too coarse to resolve what the optimum wants, so the point count is
##      doubled and the fit repeated.
##
## PARALLELISM.  N_WORKERS = 1 runs sequentially; N_WORKERS > 1 forks the walks
## with parallel::mclapply.  The walks for different K are independent, but the
## nodes WITHIN a walk are strictly sequential, since each is the warm start for
## the next.  A walk therefore cannot be split, and the speed-up saturates at
## one worker per K (four here).  Nothing draws random numbers, so forked
## workers need no seed handling.
##
## OUTPUT
##   studentt_mixtures.rds  list keyed by K, holding nu, KL and the K-column
##                     matrices w, m, v
##   studentt_mixtures.csv  one row per (K, nu): K, nu, w1..w16, m1..m16, v1..v16,
##                     blank where a mixture has fewer than 16 components
## ===========================================================================

suppressMessages({library(matrixStats); library(clue)})

## ---- user settings --------------------------------------------------------
K_SHIP     <- c(7L, 10L, 12L, 16L)          # mixtures written out
KMAX       <- max(K_SHIP)                   # seeding ladder runs K = 2..KMAX
NU_N       <- 500L                          # number of nu nodes
NU_LO      <- 2.01; NU_HI <- 500             # clamp of the Student-t nu estimator
N_WORKERS  <- 4L                             # 1 = sequential; >1 = parallel walks
GAUSS_RDS  <- "gaussian_mixtures.rds"        # optional cross-check at the largest nu

.make_nu_grid <- function(N = NU_N, lo = NU_LO, hi = NU_HI)
  rev(1/seq(1/hi, 1/lo, length.out = N))
NU_GRID    <- .make_nu_grid()
N_QUAD     <- 2000L                         # quadrature points used for fitting
N_QUAD_MAX <- 32000L                        # refinement cap
N_EVAL     <- 120000L                       # quadrature points used to evaluate KL
TOBS       <- 7775                          # only used to report nats

## Locate this file, so everything is written beside it however it was started:
## Rscript supplies --file=, while source() (the RStudio Source button) sets
## ofile on its evaluation frame.  If neither is available the working directory
## is used and no build hash is recorded; both are reported at startup.
.script_path <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(normalizePath(sub("^--file=", "", a[1])))
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of))
  }
  NA_character_
}
SCRIPT <- .script_path()
OUTDIR <- if (is.na(SCRIPT)) getwd() else dirname(SCRIPT)
HELPERS <- file.path(OUTDIR, "Table_GMKF_helpers.R")
if (!file.exists(HELPERS)) stop("cannot find ", HELPERS)

## Interval searched for the quadrature window.  A setting, not a constant:
## quad_range() lays a fixed number of points across it, so a wider interval
## coarsens that scan and shifts the window it selects, changing every value
## in the table.
QUAD_SCAN <- c(-300, 200)
source(HELPERS)

## ---- provenance -----------------------------------------------------------
## Records which version of this script produced a table.  The hash is an md5
## of this file, so any edit to it changes the hash.  The stamp is printed at
## the top of the build log and attached to the saved object, so a table found
## on disk can always be matched back to the code that produced it.  It is not
## written to the .csv, which has nowhere natural to put it.
.provenance <- function() {
  fs <- c(SCRIPT, HELPERS)
  md5 <- function(f) if (!is.na(f) && file.exists(f)) unname(tools::md5sum(f)) else NA_character_
  list(script      = if (is.na(SCRIPT)) NA_character_ else basename(SCRIPT),
       script_md5  = md5(SCRIPT),
       helpers_md5 = md5(HELPERS),
       built       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       r_version   = paste0(R.version$major, ".", R.version$minor))
}
PROV <- .provenance()


## ---- target density, centred ----------------------------------------------
## z ~ Student-t(nu) unstandardised, so z^2 ~ F(1, nu).  With x = log(z^2)
## centred by MU below, the density of x is df(u, 1, nu) * u at u = exp(x + MU).
mu_logz2 <- function(nu) digamma(0.5) - digamma(nu/2) + log(nu)
dens_target <- function(x, nu) {
  u <- exp(x + mu_logz2(nu)); d <- df(u, 1, nu) * u
  d[!is.finite(d)] <- 0; d
}
## Scale of the target: Var(x) = psi'(1/2) + psi'(nu/2) exactly.  It varies only
## from 2.5628 at nu = 2.01 to 2.2223 at nu = 500, tending to sqrt(pi^2/2) as
## nu grows, so the predictor below moves the warm start very little here.
sd_target <- function(nu) sqrt(trigamma(0.5) + trigamma(nu/2))


KL_TAIL <- 1e-16


## ---- run ------------------------------------------------------------------
## ---- output locations ------------------------------------------------------
## Console output goes to logs/, the .csv to csvTables/; both are created if
## absent.  Errors and warnings are not captured: under Rscript they go to
## stderr, in RStudio to the console pane.  Neither reaches the log file.
LOGDIR <- file.path(OUTDIR, "logs");      dir.create(LOGDIR, showWarnings = FALSE)
CSVDIR <- file.path(OUTDIR, "csvTables"); dir.create(CSVDIR, showWarnings = FALSE)
.logcon <- file(file.path(LOGDIR, "studentt_build.log"), open = "wt")
sink(.logcon, split = TRUE)
cat("Student-t mixture table for the GMKF\n")
cat(sprintf("built by %s  md5 %s\n         helpers    md5 %s\n         %s  R %s\n",
            PROV$script, PROV$script_md5, PROV$helpers_md5, PROV$built, PROV$r_version))
if (is.na(PROV$script_md5) || is.na(PROV$helpers_md5))
  warning("this script could not locate itself, so the table carries no build ",
          "hash and output goes to ", OUTDIR, call. = FALSE, immediate. = TRUE)
cat("target: centred log(z^2), z ~ Student-t(nu)\n")
cat(sprintf("nu grid: %d nodes from %.4f to %.4f, uniform in 1/nu\n",
            length(NU_GRID), min(NU_GRID), max(NU_GRID)))
cat("components:", paste(K_SHIP, collapse = ", "), "  workers:", N_WORKERS, "\n\n")

t_start <- proc.time()[["elapsed"]]
n <- length(NU_GRID); i_start <- n

## Build the starting mixture at the largest nu from scratch, by a nested
## ladder in K.  This is the only node not reached by continuation.
cat(sprintf("building the starting mixture at nu = %g (nested ladder K = 2..%d)\n",
            NU_GRID[i_start], KMAX))
start <- seed_ladder(NU_GRID[i_start])
cat(sprintf("starting mixture done, %.1f s\n\n", proc.time()[["elapsed"]] - t_start))

## Walk down the grid, one independent walk per K.
cat(sprintf("walking nu downward from %g to %g, %d walks on %d worker(s)\n\n",
            NU_GRID[i_start], NU_GRID[1], length(K_SHIP), N_WORKERS))
jobs <- lapply(K_SHIP, function(K) local({ K_ <- K; s <- start[[as.character(K)]]
  function() walk_nodes(K_, s, i_start, (i_start-1):1) }))
res <- .run_jobs(jobs)

tables <- setNames(lapply(seq_along(K_SHIP), function(j)
  .assemble(K_SHIP[j], i_start, start[[as.character(K_SHIP[j])]], res[[j]])),
  as.character(K_SHIP))
cat(sprintf("\ntotal %.1f s\n\n", proc.time()[["elapsed"]] - t_start))

## ---- write ----------------------------------------------------------------
attr(tables, "provenance") <- PROV
saveRDS(tables, file.path(OUTDIR, "studentt_mixtures.rds"))
.g17 <- function(x, nn) c(vapply(x, function(z) sprintf("%.17g", z), ""),
                          rep("", nn - length(x)))
rows <- do.call(rbind, lapply(tables, function(tb) {
  do.call(rbind, lapply(seq_along(tb$nu), function(i) {
    r <- c(tb$K, sprintf("%.10g", tb$nu[i]),
           .g17(tb$w[i, ], KMAX), .g17(tb$m[i, ], KMAX), .g17(tb$v[i, ], KMAX))
    setNames(as.data.frame(as.list(r), stringsAsFactors = FALSE),
             c("K", "nu", paste0("w", 1:KMAX), paste0("m", 1:KMAX), paste0("v", 1:KMAX)))
  })) }))
rownames(rows) <- NULL
write.csv(rows, file.path(CSVDIR, "studentt_mixtures.csv"), row.names = FALSE, quote = FALSE)

## ---- summary --------------------------------------------------------------
cat("mixtures written (centred scale)\n")
for (tb in tables)
  cat(sprintf("  K=%2d  KL from %.3e to %.3e (largest at nu=%.4f)  any non-finite: %s\n",
              tb$K, min(tb$KL), max(tb$KL), tb$nu[which.max(tb$KL)],
              any(!is.finite(c(tb$w, tb$m, tb$v)))))

## ---- cross-check at the largest nu, if the Gaussian table is available ----
## log(z^2) tends to log(chi^2_1) as nu grows, so this node should be close to
## the Gaussian mixtures without matching them exactly.  Reported only.
gp <- file.path(OUTDIR, GAUSS_RDS)
if (file.exists(gp)) {
  gz <- readRDS(gp); iN <- length(NU_GRID)
  cat(sprintf("\ncross-check at nu = %g against %s (a tolerance, not an identity)\n",
              NU_GRID[iN], GAUSS_RDS))
  for (Kc in names(tables)) {
    if (is.null(gz[[Kc]])) { cat(sprintf("  K=%2s  not present in %s\n", Kc, GAUSS_RDS)); next }
    tb <- tables[[Kc]]; z <- gz[[Kc]]
    p <- track(list(w = z$w, m = z$m, v = z$v),
               list(w = tb$w[iN, ], m = tb$m[iN, ], v = tb$v[iN, ]))
    cat(sprintf("  K=%2s  KL ratio %.6f   max rel dw %.2e   max |dm| %.2e   max rel dv %.2e\n",
                Kc, tb$KL[iN]/z$KL, max(abs(p$w/z$w - 1)),
                max(abs(p$m - z$m)), max(abs(p$v/z$v - 1))))
  }
} else {
  cat(sprintf("\ncross-check at nu = %g skipped: %s not found\n", max(NU_GRID), GAUSS_RDS))
}
cat("\nwritten: studentt_mixtures.rds  studentt_mixtures.csv\n")

sink(); close(.logcon)
