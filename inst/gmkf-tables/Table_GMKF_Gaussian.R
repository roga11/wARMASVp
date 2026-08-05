## ===========================================================================
## Table_GMKF_Gaussian.R
##
## Builds KL-optimal Gaussian-mixture approximations to the density of the
## centred log squared innovation
##
##     x = log(z_t^2) - E[log(z_t^2)],     z_t ~ N(0,1),
##
## for K = 7, 10, 12 and 16 components.  This is the measurement-error density
## of the linearised SV(p) observation equation
##
##     y*_t = log(y_t^2) - mu = h_t + x_t,   mu = log(sigma_y^2) + E[log(z_t^2)],
##
## and the mixture is what the Gaussian-mixture Kalman filter (GMKF) uses in
## place of it.
##
## CENTRING.  Everything here is on the CENTRED scale, matching y*_t above and
## matching the mixtures the package fits for Student-t and GED.  For the
## Gaussian case E[log(z_t^2)] = psi(1/2) + log 2 = -1.2703628455, so the
## fitted component means have weighted mean 0 by construction, not -1.27.
## (Kim, Shephard and Chib (1998) Table 4 is published on a centred scale too,
## but with the offset rounded to 1.2704.)
##
## NO nu GRID.  The target is parameter-free -- log(chi^2_1) has no shape
## parameter -- so the deliverable is four mixtures, not a table indexed by nu
## and interpolated.  Nothing here interpolates.
##
## TWO GRIDS, and they are different things:
##   QUADRATURE POINTS  the values of x at which the target density is
##                      evaluated to form the objective.  2000 of them here.
##   NODES              tabulated values of a shape parameter nu.  The
##                      Gaussian target has none.
##
## METHOD
##   1. Deterministic quadrature.  The target is evaluated at 2000 quadrature
##      points spanning all but ~1e-13 of its mass, giving weights p_i, and the
##      objective is
##          F(w,m,v) = - sum_i p_i log g(x_i),        F = KL + constant,
##      minimised over the mixture g.  The mixture is fitted to a DENSITY, not
##      to data: the quadrature points and their weights play the role of a
##      sample.  There is no sampling and no RNG anywhere, so the result is
##      bit-reproducible and carries no seed.  This is a one-off table-building
##      computation; nothing at filter runtime repeats it.
##   2. Unconstrained gauge-fixed chart.  w = softmax(eta) with eta_1 = 0 held
##      fixed, v = exp(ell).  Parameters th = (eta_2..eta_K, m_1..m_K,
##      ell_1..ell_K), dimension 3K-1.  Gradient and Hessian are analytic.
##   3. Weighted EM, then a trust-region Newton polish.  EM is used only as a
##      globaliser: it is monotone and cannot diverge, so it reliably carries a
##      crude start into a good basin, but its convergence rate approaches 1
##      and it stalls long before reaching the optimum.  The trust-region step
##      then converges it, using exact eigendecomposition subproblem solves so
##      that indefinite Hessians and the hard case are handled rather than
##      damped away.
##   4. Nested continuation in K.  The ladder runs K = 2..16; each K is started
##      by splitting each component of the polished K-1 fit in turn and keeping
##      whichever candidate polishes to the lowest objective.  Splitting at
##      m +/- eps*sd with weights w/2 and variance v*(1-eps^2) preserves the
##      second moment, so the start is density-comparable to the K-1 solution.
##      Only K = 7, 10, 12, 16 are written out; the intermediate rungs exist to
##      seed them.
##   5. Adaptive quadrature refinement.  A variance floor vmin = (2*Delta)^2 is
##      imposed, Delta being the quadrature spacing.  It is needed because the
##      DISCRETISED objective is unbounded: a component whose variance shrinks
##      to zero while centred on a quadrature point drives F to -infinity,
##      which is an artefact of the discretisation and not a real optimum.  If
##      a fit lands ON that floor the grid is too coarse to resolve what the
##      optimum wants, so the point count is doubled and the fit repeated --
##      the component is legitimate and the grid is at fault.  The Gaussian
##      target never approaches the floor (the tightest component sits several
##      multiples above it), so this never triggers here.
##
## OUTPUT
##   gaussian_mixtures.rds  list keyed by K, each with w, m, v and KL
##   gaussian_mixtures.csv  one row per K: K, w1..w16, m1..m16, v1..v16
##                          (blank where a mixture has fewer than 16 components)
## ===========================================================================

suppressMessages(library(matrixStats))   # rowMaxs

## ---- user settings --------------------------------------------------------
K_SHIP      <- c(7L, 10L, 12L, 16L)  # mixtures written out
KMAX        <- max(K_SHIP)           # ladder runs K = 2..KMAX
N_QUAD      <- 2000L                 # quadrature points used for fitting
N_QUAD_MAX  <- 32000L                # refinement cap
N_EVAL      <- 120000L               # quadrature points used to evaluate KL
TOBS        <- 7775                  # only used to report nats

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
QUAD_SCAN <- c(-200, 120)
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
## X = log(u), u ~ chi^2_1, has density dchisq(e^X, 1) * e^X and mean MU_LOGZ2.
## The centred variable is x = X - MU_LOGZ2, so f(x) = f_X(x + MU_LOGZ2).
MU_LOGZ2   <- digamma(0.5) + log(2)          # -1.2703628455 = E[log(z^2)]
## The nu argument is unused -- log(chi^2_1) has no shape parameter -- and is
## present only so the shared quadrature code takes the same form for every
## error family.
dens_target <- function(x, nu = NULL) { u <- exp(x + MU_LOGZ2); dchisq(u, 1) * u }


KL_TAIL <- 1e-16


## ---- the ladder: K = 2 .. KMAX, nested ------------------------------------
build_ladder <- function(Kmax = KMAX, n_quad = N_QUAD, n_quad_max = N_QUAD_MAX) {
  qd <- make_quad(NULL, n_quad); prev <- NULL; out <- list()
  for (K in 2:Kmax) {
    t0 <- proc.time()[["elapsed"]]
    st <- fit_step(K, qd, prev)
    ## a fit sitting ON the variance floor means the quadrature is too coarse
    nq <- length(qd$x)
    while (min(sqrt(st$v)) <= 2 * qd$delta * 1.001 && nq < n_quad_max) {
      nq <- 2L * nq; qd <- make_quad(NULL, nq)
      st <- fit_step(K, qd, prev)
    }
    KL <- kl_eval(NULL, st$w, st$m, st$v)
    out[[as.character(K)]] <- list(K = K, w = st$w, m = st$m, v = st$v,
                                   F = st$F, KL = KL, nats = TOBS * KL,
                                   n_quad = length(qd$x), which = st$which,
                                   gnorm = st$gnorm, stop_reason = st$stop_reason,
                                   n_restart = st$n_restart,
                                   secs = proc.time()[["elapsed"]] - t0)
    cat(sprintf("K=%2d  KL=%.6e  nats=%.6f  quad=%5d  %-8s  |g|=%.2e  %-16s  %5.1fs\n",
                K, KL, TOBS * KL, length(qd$x), st$which, st$gnorm,
                st$stop_reason, out[[as.character(K)]]$secs))
    flush.console()
    prev <- list(w = st$w, m = st$m, v = st$v, F = st$F)
  }
  out
}

## ---- run ------------------------------------------------------------------
## ---- output locations ------------------------------------------------------
## Console output goes to logs/, the .csv to csvTables/; both are created if
## absent.  Errors and warnings are not captured: under Rscript they go to
## stderr, in RStudio to the console pane.  Neither reaches the log file.
LOGDIR <- file.path(OUTDIR, "logs");      dir.create(LOGDIR, showWarnings = FALSE)
CSVDIR <- file.path(OUTDIR, "csvTables"); dir.create(CSVDIR, showWarnings = FALSE)
.logcon <- file(file.path(LOGDIR, "run_gaussian.log"), open = "wt")
sink(.logcon, split = TRUE)
cat(sprintf("built by %s  md5 %s\n         helpers    md5 %s\n         %s  R %s\n",
            PROV$script, PROV$script_md5, PROV$helpers_md5, PROV$built, PROV$r_version))
if (is.na(PROV$script_md5) || is.na(PROV$helpers_md5))
  warning("this script could not locate itself, so the table carries no build ",
          "hash and output goes to ", OUTDIR, call. = FALSE, immediate. = TRUE)
cat("Gaussian GMKF mixtures -- centred log(z^2), E[log(z^2)] = ",
    sprintf("%.10f", MU_LOGZ2), "\n", sep = "")
cat("nested ladder K = 2 ..", KMAX, "; written out:", paste(K_SHIP, collapse = ", "), "\n\n")

t_start <- proc.time()[["elapsed"]]
ladder  <- build_ladder()
cat(sprintf("\ntotal %.1f s\n\n", proc.time()[["elapsed"]] - t_start))

tables <- ladder[as.character(K_SHIP)]

## ---- write ----------------------------------------------------------------
attr(tables, "provenance") <- PROV
saveRDS(tables, file.path(OUTDIR, "gaussian_mixtures.rds"))

## One row per K: K, then the K weights, means and variances, padded with blanks
## out to KMAX columns.  17 significant digits; the .rds holds the same values
## as exact doubles.
.g17 <- function(x, n) c(vapply(x, function(z) sprintf("%.17g", z), ""),
                         rep("", n - length(x)))
rows <- do.call(rbind, lapply(tables, function(z) {
  r <- c(z$K, .g17(z$w, KMAX), .g17(z$m, KMAX), .g17(z$v, KMAX))
  setNames(as.data.frame(as.list(r), stringsAsFactors = FALSE),
           c("K", paste0("w", 1:KMAX), paste0("m", 1:KMAX), paste0("v", 1:KMAX)))
}))
rownames(rows) <- NULL
write.csv(rows, file.path(CSVDIR, "gaussian_mixtures.csv"),
          row.names = FALSE, quote = FALSE)

## ---- summary --------------------------------------------------------------
cat("written mixtures (centred scale)\n")
for (z in tables)
  cat(sprintf("  K=%2d  KL=%.6e  nats(T=%d)=%.6f  sum(w*m)=%+.3e  min(w)=%.3e  min(sd)=%.4f\n",
              z$K, z$KL, TOBS, z$nats, sum(z$w * z$m), min(z$w), min(sqrt(z$v))))
cat("\nwritten:\n  ", file.path(OUTDIR, "gaussian_mixtures.rds"),
    "\n  ", file.path(CSVDIR, "gaussian_mixtures.csv"), "\n", sep = "")

sink(); close(.logcon)
