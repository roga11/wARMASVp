## ===========================================================================
## Table_GMKF_GED.R
##
## Builds KL-optimal Gaussian-mixture approximations to the density of the
## centred log squared innovation
##
##     x = log(z_t^2) - E[log(z_t^2)],     z_t ~ standardised GED(nu),
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
## mixtures the package fits for the other error families.  For standardised
## GED(nu),  E[log(z_t^2)] = (2/nu) psi(1/nu) + lgamma(1/nu) - lgamma(3/nu),
## which at nu = 2 equals psi(1/2) + log 2 = -1.2703628455, as it must, since
## GED(2) is the standard normal.
##
## NU GRID.  nu in [0.1, 20], 200 nodes, spaced uniformly in log(nu) within each
## of [0.1, 2] and [2, 20] -- about 2.7% per node.  The range is set by the estimator, which clamps
## GED nu to exactly that interval (wARMASVp/src/rootfinding.cpp: fixed_lower =
## 0.1, fixed_upper = 20.0; the boundary warnings in R/svp_internal.R fire at
## nu <= 0.1 and nu >= 20).  A table over a narrower range would leave reachable
## estimates untabulated.
##
## Why log(nu) and not nu.  How fast the target moves per unit nu is wildly
## uneven: the Hellinger distance between adjacent targets on a uniform-nu grid
## has a worst/median ratio of 131.7, with the worst step at the very bottom
## (nu 0.1 -> 0.2).  Half the family's total path length lies below nu = 0.88,
## which a uniform-nu grid barely resolves.  Spacing uniformly in log(nu) brings
## that ratio to 1.1 -- within 18% of the arclength-optimal grid -- and cuts the
## worst step 25-fold.  Continuation is only as good as its step size, so this
## is what makes each warm start a genuinely small perturbation.
## The grid is DEFINED in log(nu) and exponentiated; everything stored, printed
## and returned is plain nu.
## nu = 2 is forced onto the grid, because GED(2) is exactly the standard normal
## and that node can therefore be checked against an independently built table.
## A single uniform log grid on [0.1, 20] cannot contain nu = 2 at any N, so the
## grid is split at 2 and the node counts chosen to match the two step sizes as
## closely as an integer split allows (0.98% apart at N = 200, 0.05% at N = 2000).
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
##      every component carries real weight.  At the small-nu end the
##      approximation is easy and some components become redundant, which
##      leaves the fit only loosely pinned down -- a poor place to begin a
##      sequence of nodes that every later node depends on.
##   4b. Optional check at nu = 2.  A standardised GED(2) variable IS standard
##      normal, so log(z^2) at nu = 2 is exactly log(chi^2_1) -- the density
##      Table_GMKF_Gaussian.R solves by a separate route.  If its output is
##      present, the walk's nu = 2 node is compared against it.  This is an
##      independent check on the whole walk, since nu = 2 is reached after
##      roughly ninety continuation steps.  It is skipped if that file is
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
##      doubled and the fit repeated.  GED at large nu is the case that needs
##      this: the density develops a hard edge as nu grows.
##
## PARALLELISM.  N_WORKERS = 1 runs sequentially; N_WORKERS > 1 forks the walks
## with parallel::mclapply.  The walks for different K are independent, but the
## nodes WITHIN a walk are strictly sequential, since each is the warm start for
## the next.  A walk therefore cannot be split, and the speed-up saturates at
## one worker per K (four here).  Nothing draws random numbers, so forked
## workers need no seed handling.
##
## OUTPUT
##   ged_mixtures.rds  list keyed by K, each holding nu, KL and the K-column
##                     matrices w, m, v
##   ged_mixtures.csv  one row per (K, nu): K, nu, w1..w16, m1..m16, v1..v16,
##                     blank where a mixture has fewer than 16 components
## ===========================================================================

suppressMessages({library(matrixStats); library(clue)})

## ---- user settings --------------------------------------------------------
K_SHIP     <- c(7L, 10L, 12L, 16L)          # mixtures written out
KMAX       <- max(K_SHIP)                   # seeding ladder runs K = 2..KMAX
NU_N       <- 200L                          # number of nu nodes
NU_LO      <- 0.1; NU_HI <- 20               # clamp of the GED nu estimator
NU_ANCHOR  <- 2                              # forced onto the grid: GED(2) == N(0,1)
N_WORKERS  <- 4L                             # 1 = sequential; >1 = parallel walks
GAUSS_RDS  <- "gaussian_mixtures.rds"        # optional cross-check at nu = 2

## Grid: uniform in log(nu) within each of [NU_LO, NU_ANCHOR] and
## [NU_ANCHOR, NU_HI], with the node counts split so the two step sizes match
## as closely as an integer split allows.  This puts BOTH endpoints and
## nu = 2 exactly on the grid while keeping the spacing effectively uniform in
## log(nu): at N = 200 the two steps differ by 0.98%, at N = 2000 by 0.05%.
## A single uniform log grid on [0.1, 20] cannot contain nu = 2 for any N,
## because that needs log(20)/log(200) = 0.56541 times (N-1) to be an integer.
.make_nu_grid <- function(N = NU_N, lo = NU_LO, hi = NU_HI, anchor = NU_ANCHOR) {
  a <- log(lo); b <- log(anchor); c_ <- log(hi)
  n1 <- round((N - 1L) * (b - a) / (c_ - a)); n2 <- (N - 1L) - n1
  c(exp(seq(a, b, length.out = n1 + 1L)), exp(seq(b, c_, length.out = n2 + 1L))[-1])
}
NU_GRID    <- .make_nu_grid()
IDX_ANCHOR <- which.min(abs(NU_GRID - NU_ANCHOR))
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
## z standardised GED(nu):  f_z(z) = (nu / (2 a Gamma(1/nu))) exp(-|z/a|^nu),
## a = sqrt(Gamma(1/nu)/Gamma(3/nu)) giving Var(z) = 1.  With u = z^2 and
## x = log(u) - MU(nu), the density of x is f_z(sqrt(u)) * u / sqrt(u).
mu_logz2 <- function(nu) (2/nu)*digamma(1/nu) + lgamma(1/nu) - lgamma(3/nu)
## Scale of the target.  Writing G ~ Gamma(1/nu, 1), the centred variable is
## exactly  x = (2/nu) * (log G - psi(1/nu)),  so sd(x) = (2/nu)*sqrt(psi'(1/nu)).
## The whole nu-dependence is a KNOWN scale factor times a log-Gamma shape.
sd_target <- function(nu) (2/nu)*sqrt(trigamma(1/nu))
dens_target <- function(x, nu) {
  a  <- exp(0.5*(lgamma(1/nu) - lgamma(3/nu)))
  u2 <- exp(x + mu_logz2(nu)); u <- sqrt(u2)
  cg <- nu/(2*a*gamma(1/nu))
  2*cg*exp(-(u/a)^nu) * (1/(2*u)) * u2
}


KL_TAIL <- 1e-16


## ---- run ------------------------------------------------------------------
## ---- output locations ------------------------------------------------------
## Console output goes to logs/, the .csv to csvTables/; both are created if
## absent.  Errors and warnings are not captured: under Rscript they go to
## stderr, in RStudio to the console pane.  Neither reaches the log file.
LOGDIR <- file.path(OUTDIR, "logs");      dir.create(LOGDIR, showWarnings = FALSE)
CSVDIR <- file.path(OUTDIR, "csvTables"); dir.create(CSVDIR, showWarnings = FALSE)
.logcon <- file(file.path(LOGDIR, "ged_build.log"), open = "wt")
sink(.logcon, split = TRUE)
cat("GED mixtures for the Gaussian-mixture Kalman filter\n")
cat(sprintf("built by %s  md5 %s\n         helpers    md5 %s\n         %s  R %s\n",
            PROV$script, PROV$script_md5, PROV$helpers_md5, PROV$built, PROV$r_version))
if (is.na(PROV$script_md5) || is.na(PROV$helpers_md5))
  warning("this script could not locate itself, so the table carries no build ",
          "hash and output goes to ", OUTDIR, call. = FALSE, immediate. = TRUE)
cat("target: centred log(z^2), z standardised GED(nu)\n")
cat(sprintf("nu grid: %d nodes from %.4f to %.4f, uniform in log(nu) each side of %g\n",
            length(NU_GRID), min(NU_GRID), max(NU_GRID), NU_ANCHOR))
cat(sprintf("  %.4f%% per node below nu=%g, %.4f%% above\n",
            100*(NU_GRID[2]/NU_GRID[1]-1), NU_ANCHOR,
            100*(NU_GRID[length(NU_GRID)]/NU_GRID[length(NU_GRID)-1]-1)))
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
saveRDS(tables, file.path(OUTDIR, "ged_mixtures.rds"))
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
write.csv(rows, file.path(CSVDIR, "ged_mixtures.csv"), row.names = FALSE, quote = FALSE)

## ---- summary --------------------------------------------------------------
cat("mixtures written (centred scale)\n")
for (tb in tables)
  cat(sprintf("  K=%2d  KL from %.3e to %.3e (largest at nu=%.4f)  any non-finite: %s\n",
              tb$K, min(tb$KL), max(tb$KL), tb$nu[which.max(tb$KL)],
              any(!is.finite(c(tb$w, tb$m, tb$v)))))

## ---- cross-check at nu = 2, if the Gaussian table is available -------------
## A standardised GED(2) variable is standard normal, so this node solves the
## same density as Table_GMKF_Gaussian.R by an independent route.  Reported
## only; nothing above depends on it.
gp <- file.path(OUTDIR, GAUSS_RDS)
if (file.exists(gp)) {
  gz <- readRDS(gp); i2 <- which.min(abs(NU_GRID - NU_ANCHOR))
  cat(sprintf("\ncross-check at nu = %g against %s\n", NU_GRID[i2], GAUSS_RDS))
  for (Kc in names(tables)) {
    if (is.null(gz[[Kc]])) { cat(sprintf("  K=%2s  not present in %s\n", Kc, GAUSS_RDS)); next }
    tb <- tables[[Kc]]; z <- gz[[Kc]]
    p <- track(list(w = z$w, m = z$m, v = z$v),
               list(w = tb$w[i2, ], m = tb$m[i2, ], v = tb$v[i2, ]))
    cat(sprintf("  K=%2s  KL ratio %.8f   max rel dw %.2e   max |dm| %.2e   max rel dv %.2e\n",
                Kc, tb$KL[i2]/z$KL, max(abs(p$w/z$w - 1)),
                max(abs(p$m - z$m)), max(abs(p$v/z$v - 1))))
  }
} else {
  cat(sprintf("\ncross-check at nu = %g skipped: %s not found\n", NU_ANCHOR, GAUSS_RDS))
}
cat("\nwritten: ged_mixtures.rds  ged_mixtures.csv\n")

sink(); close(.logcon)
