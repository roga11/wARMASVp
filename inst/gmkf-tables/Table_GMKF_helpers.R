## ===========================================================================
## Table_GMKF_helpers.R
##
## Shared machinery for the three table builders -- Table_GMKF_Gaussian.R,
## Table_GMKF_GED.R and Table_GMKF_StudentT.R -- and for the diagnostics in
## diagnostics/Certify_GMKF_Tables.R.  Sourcing this file defines functions and
## nothing else: it starts no build, writes no file and creates no directory.
##
## Each build script supplies what actually differs between the error families
## and then sources this file:
##
##   dens_target(x, nu)   the density of the centred log squared innovation.
##                        Written with a nu argument even where there is no
##                        shape parameter, so the quadrature code below is
##                        common to all three; the Gaussian version ignores it.
##   QUAD_SCAN            the interval searched for the quadrature window.  It
##                        is a setting, not a constant, because quad_range()
##                        lays a fixed number of points across it -- widening
##                        the interval coarsens that scan and moves the window
##                        it selects, which would change every fitted value.
##   NU_GRID              the tabulated shape parameters, if the family has one
##   sd_target(nu)        scale of the target, for the continuation predictor
##   N_QUAD, N_QUAD_MAX, N_EVAL, KL_TAIL, TOBS, K_SHIP, KMAX, N_WORKERS
##
## The mathematics is documented in the build scripts; this file holds the
## implementation only.
## ===========================================================================

suppressMessages({library(matrixStats); library(clue)})

## ---- compensated summation ------------------------------------------------
.neumaier <- function(z) {
  s <- 0.0; cmp <- 0.0
  for (i in seq_along(z)) {
    zi <- z[i]; t <- s + zi
    if (abs(s) >= abs(zi)) cmp <- cmp + ((s - t) + zi) else cmp <- cmp + ((zi - t) + s)
    s <- t
  }
  s + cmp
}

.csum <- function(z, block = 1024L) {
  n <- length(z)
  if (n == 0L) return(0.0)
  if (anyNA(z) || any(!is.finite(z))) return(sum(z))
  if (n <= block) return(.neumaier(z))
  nb <- as.integer(ceiling(n / block))
  length(z) <- nb * block
  z[is.na(z)] <- 0
  .csum(colSums(matrix(z, nrow = block)), block)
}

## ---- quadrature -----------------------------------------------------------
## The FITTING window leaves ~1e-13 of the mass outside, and is kept tight on
## purpose: with a fixed point count a wider window means a larger Delta, which
## raises the variance floor.  The EVALUATION window is wider (1e-16) because
## it has its own point count and nothing is coarsened by widening it.
quad_range <- function(nu = NULL, lo = QUAD_SCAN[1], hi = QUAD_SCAN[2],
                       n = 200000L, tail = 1e-13) {
  x <- seq(lo, hi, length.out = n); f <- dens_target(x, nu); f[!is.finite(f)] <- 0
  cf <- cumsum(f); cf <- cf / cf[length(cf)]
  c(x[max(1L, which(cf > tail)[1] - 1L)], x[min(n, which(cf > 1 - tail)[1] + 1L)])
}

make_quad <- function(nu = NULL, n_quad) {
  rg <- quad_range(nu)
  x  <- seq(rg[1], rg[2], length.out = n_quad)
  p  <- dens_target(x, nu); p[!is.finite(p)] <- 0; p <- p / sum(p)
  list(x = x, p = p, delta = x[2] - x[1], range = rg)
}

## ---- chart ----------------------------------------------------------------
pack   <- function(w, m, v) c(log(w[-1] / w[1]), m, log(v))

unpack <- function(th, K) {
  eta <- c(0, th[1:(K - 1)]); w <- exp(eta - max(eta)); w <- w / sum(w)
  list(w = w, m = th[K:(2 * K - 1)], v = exp(th[(2 * K):(3 * K - 1)]))
}

## ---- mixture log-density, sufficient statistics, objective, gradient ------
log_mix <- function(x, w, m, v) {
  K <- length(w); L <- matrix(0, length(x), K)
  for (k in 1:K) L[, k] <- log(w[k]) - 0.5 * log(2 * pi * v[k]) - (x - m[k])^2 / (2 * v[k])
  mx <- rowMaxs(L); mx + log(rowSums(exp(L - mx)))
}

suffstats <- function(qd, w, m, v) {
  x <- qd$x; p <- qd$p; K <- length(w); L <- matrix(0, length(x), K)
  for (k in 1:K) L[, k] <- log(w[k]) - 0.5 * log(2 * pi * v[k]) - (x - m[k])^2 / (2 * v[k])
  mx <- rowMaxs(L); lg <- mx + log(rowSums(exp(L - mx))); R <- exp(L - lg)
  S <- numeric(K); for (k in 1:K) S[k] <- sum(p * R[, k] * (x - m[k])^2)
  list(W = colSums(p * R), A = colSums(p * R * x), S = S, loglik = sum(p * lg))
}

objF <- function(qd, w, m, v) -sum(qd$p * log_mix(qd$x, w, m, v))

grad_chart <- function(qd, w, m, v) {
  ss <- suffstats(qd, w, m, v)
  c((w - ss$W)[-1], -(ss$A - m * ss$W) / v, -(ss$S - v * ss$W) / (2 * v))
}

Fth <- function(qd, th, K) { q <- unpack(th, K); objF(qd, q$w, q$m, q$v) }

Gth <- function(qd, th, K) { q <- unpack(th, K); grad_chart(qd, q$w, q$m, q$v) }

## ---- analytic Hessian in the chart ---------------------------------------
hess_chart <- function(qd, w, m, v) {
  x <- qd$x; p <- qd$p; K <- length(w)
  L <- matrix(0, length(x), K)
  for (k in 1:K) L[, k] <- log(w[k]) - 0.5 * log(2 * pi * v[k]) - (x - m[k])^2 / (2 * v[k])
  mx <- rowMaxs(L); lg <- mx + log(rowSums(exp(L - mx))); R <- exp(L - lg)
  B <- Q <- matrix(0, length(x), K)
  for (k in 1:K) { B[, k] <- (x - m[k]) / v[k]; Q[, k] <- (x - m[k])^2 / (2 * v[k]) - 0.5 }
  Teta <- sweep(R, 2, w, "-")[, -1, drop = FALSE]
  Tall <- cbind(Teta, R * B, R * Q)
  H <- crossprod(Tall * sqrt(p), Tall * sqrt(p))
  ie <- seq_len(K - 1); im <- K:(2 * K - 1); il <- (2 * K):(3 * K - 1)
  for (a in 1:(K - 1)) { ka <- a + 1
    for (b in 1:(K - 1)) { kb <- b + 1
      cab <- (if (ka == kb) (R[, ka] - w[ka]) else 0) - w[ka] * R[, kb] -
             w[kb] * R[, ka] + 2 * w[ka] * w[kb]
      H[ie[a], ie[b]] <- H[ie[a], ie[b]] - sum(p * cab) }
    for (l in 1:K) {
      cml <- R[, l] * B[, l] * ((ka == l) - w[ka])
      cll <- R[, l] * Q[, l] * ((ka == l) - w[ka])
      H[ie[a], im[l]] <- H[ie[a], im[l]] - sum(p * cml); H[im[l], ie[a]] <- H[ie[a], im[l]]
      H[ie[a], il[l]] <- H[ie[a], il[l]] - sum(p * cll); H[il[l], ie[a]] <- H[ie[a], il[l]] } }
  for (k in 1:K) {
    H[im[k], im[k]] <- H[im[k], im[k]] - sum(p * R[, k] * (B[, k]^2 - 1 / v[k]))
    H[im[k], il[k]] <- H[im[k], il[k]] - sum(p * R[, k] * B[, k] * (Q[, k] - 1))
    H[il[k], im[k]] <- H[im[k], il[k]]
    H[il[k], il[k]] <- H[il[k], il[k]] - sum(p * R[, k] * (Q[, k]^2 - Q[, k] - 0.5)) }
  (H + t(H)) / 2
}

## ---- weighted EM (globaliser; deterministic, no sampling) -----------------
em_run <- function(qd, w, m, v, max_iter = 2000L, vmin = 0) {
  for (it in 1:max_iter) {
    ss <- suffstats(qd, w, m, v)
    if (!all(is.finite(ss$W)) || !all(is.finite(ss$A)) || !all(is.finite(ss$S))) break
    W <- pmax(ss$W, 1e-300)
    w <- W / sum(W); m <- ss$A / W
    v <- pmax(ss$S / W, vmin)
    if (!all(is.finite(w)) || !all(is.finite(m)) || !all(is.finite(v))) break
  }
  list(w = w, m = m, v = v)
}

## ---- trust-region subproblem, solved exactly ------------------------------
tr_solve <- function(g, H, radius) {
  eg <- eigen(H, symmetric = TRUE); lam1 <- min(eg$values)
  Q <- eg$vectors; ghat <- as.numeric(crossprod(Q, g))
  scale <- max(abs(eg$values), .Machine$double.eps)
  dv <- eg$values - lam1
  pnorm_u <- function(u) sqrt(sum((ghat / (dv + u))^2))
  pvec_u  <- function(u) -as.numeric(Q %*% (ghat / (dv + u)))
  mval    <- function(p) sum(g * p) + 0.5 * as.numeric(crossprod(p, H %*% p))
  if (lam1 > 0) {
    p <- -as.numeric(Q %*% (ghat / eg$values))
    if (sqrt(sum(p^2)) <= radius) return(list(p = p, bdry = FALSE))
  }
  lo <- max(0, -lam1) + 1e-14
  if (pnorm_u(lam1 + lo) < radius) {          # hard case
    idx  <- which(eg$values <= lam1 + 1e-12 * scale)
    keep <- abs(eg$values + lo) > 1e-12 * scale
    p0 <- -as.numeric(Q[, keep, drop = FALSE] %*% (ghat[keep] / (eg$values[keep] + lo)))
    u  <- Q[, idx[1]]
    aa <- sum(u^2); bb <- 2 * sum(p0 * u); cc <- sum(p0^2) - radius^2
    sq <- sqrt(max(0, bb^2 - 4 * aa * cc))
    p_p <- p0 + ((-bb + sq) / (2 * aa)) * u
    p_m <- p0 + ((-bb - sq) / (2 * aa)) * u
    return(list(p = if (mval(p_m) < mval(p_p)) p_m else p_p, bdry = TRUE))
  }
  u_lo <- max(0, lam1) + 1e-300
  u_hi <- u_lo + max(1, sqrt(sum(g^2)) / radius)
  while (pnorm_u(u_hi) > radius && u_hi < 1e300) u_hi <- u_hi * 2
  for (i in 1:200) {
    u_mid <- if (u_lo <= 0) u_hi / 2 else sqrt(u_lo * u_hi)
    if (pnorm_u(u_mid) > radius) u_lo <- u_mid else u_hi <- u_mid
    if (u_hi - u_lo <= 1e-15 * u_hi) break
    if (u_lo > 0 && u_hi / u_lo - 1 < 1e-14) break
  }
  u <- 0.5 * (u_lo + u_hi); p <- pvec_u(u)
  if (!is.finite(sum(p^2)) || sqrt(sum(p^2)) > radius) { u <- u_hi; p <- pvec_u(u) }
  pn <- sqrt(sum(p^2))
  if (is.finite(pn) && pn > radius) p <- p * (radius / pn)
  list(p = p, bdry = TRUE)
}

## ---- trust-region polish --------------------------------------------------
newton_tr <- function(qd, w, m, v, max_it = 1000L, vmin = 0,
                      radius0 = 1, radius_max = 1e3) {
  K <- length(w); th <- pack(w, m, v)
  f0 <- Fth(qd, th, K); g0 <- Gth(qd, th, K); rad <- radius0
  stop_reason <- "max_it"
  for (it in 1:max_it) {
    qq <- unpack(th, K); H <- hess_chart(qd, qq$w, qq$m, qq$v)
    eg <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (sqrt(sum(g0^2)) < 1e-14 && min(eg) >= -1e-10 * max(abs(eg))) {
      stop_reason <- "second_order"; break }
    sol <- tr_solve(g0, H, rad); p <- sol$p
    pred <- -(sum(g0 * p) + 0.5 * as.numeric(crossprod(p, H %*% p)))
    if (is.finite(pred) && pred <= 8 * .Machine$double.eps * abs(f0) &&
        min(eg) >= -1e-10 * max(abs(eg))) { stop_reason <- "resolution"; break }
    tn <- th + p
    if (vmin > 0) tn[(2 * K):(3 * K - 1)] <- pmax(tn[(2 * K):(3 * K - 1)], log(vmin))
    fn <- tryCatch(Fth(qd, tn, K), error = function(e) NA_real_)
    ared <- f0 - fn
    rho <- if (!is.finite(fn) || pred <= 0) -Inf else ared / pred
    if (rho > 0.10 && is.finite(fn) && ared > 0) {
      th <- tn; f0 <- fn; g0 <- Gth(qd, th, K)
      if (rho > 0.75 && sol$bdry) rad <- min(2 * rad, radius_max)
    } else {
      rad <- min(rad, sqrt(sum(p^2))) / 4
    }
    if (rad < 1e-14) { stop_reason <- "radius_collapsed"; break }
  }
  qq <- unpack(th, K)
  list(w = qq$w, m = qq$m, v = qq$v, F = f0,
       gnorm = sqrt(sum(Gth(qd, th, K)^2)), stop_reason = stop_reason)
}

## ---- one fit: cold, split-from-(K-1), or warm-from-neighbouring-nu --------
fit_step <- function(K, qd, prev = NULL, em_budget = 2000L, max_it = 1000L,
                     n_restart = 4L) {
  vmin <- (2 * qd$delta)^2

  polish <- function(w, m, v) {
    e <- em_run(qd, w, m, v, max_iter = em_budget, vmin = vmin)
    newton_tr(qd, e$w, e$m, e$v, vmin = vmin, max_it = max_it)
  }
  ## max_it and radius_collapsed are EXITS, not convergence, and either can
  ## stop away from a stationary point.  Re-enter warm with a fresh radius.
  converge <- function(p) {
    nres <- 0L
    while (nres < n_restart &&
           isTRUE(p$stop_reason %in% c("max_it", "radius_collapsed"))) {
      p2 <- newton_tr(qd, p$w, p$m, p$v, vmin = vmin, max_it = max_it)
      if (!is.finite(p2$F) || p2$F > p$F + 8 * .Machine$double.eps * abs(p$F)) break
      nres <- nres + 1L; done <- isTRUE(p2$stop_reason %in% c("resolution", "second_order"))
      p <- p2
      if (done) break
    }
    p$n_restart <- nres; p
  }

  mode <- if (is.null(prev)) "cold" else
          if (length(prev$w) == K - 1L) "split_K" else "warm_nu"

  if (mode == "cold") {
    br  <- c(0, cumsum(qd$p))
    idx <- sapply((seq_len(K) - 0.5) / K, function(q) which.min(abs(br - q)))
    b <- polish(rep(1 / K, K), qd$x[idx],
                rep(max(vmin, diff(range(qd$x))^2 / (12 * K^2) + 1), K))
  } else if (mode == "warm_nu") {
    b <- polish(prev$w, prev$m, pmax(prev$v, vmin))
  } else {
    ## Nested split.  Children at m +/- eps*sd with weights w/2 and variance
    ## v*(1-eps^2) have second moment v about m, so the split preserves the
    ## second moment.  Candidates are ranked by the POLISHED objective; ranking
    ## by a short EM run instead favours whichever split perturbs the fit
    ## least, which is the opposite of the one that pays off after polishing.
    eps <- 0.25
    cands <- lapply(seq_along(prev$w), function(j) {
      sd_ <- sqrt(prev$v[j])
      polish(c(prev$w[-j], prev$w[j] / 2, prev$w[j] / 2),
             c(prev$m[-j], prev$m[j] - eps * sd_, prev$m[j] + eps * sd_),
             c(prev$v[-j], rep(max(vmin, prev$v[j] * (1 - eps^2)), 2)))
    })
    b <- cands[[which.min(sapply(cands, function(z) z$F))]]
  }
  b <- converge(b)
  ## order by mean only when there is no predecessor to match against; the
  ## caller re-orders by optimal assignment when walking nu.
  if (mode != "warm_nu") { o <- order(b$m); b$w <- b$w[o]; b$m <- b$m[o]; b$v <- b$v[o] }
  b$mode <- mode
  b
}

## ---- optimal assignment between adjacent nu -------------------------------
track <- function(prev, cur) {
  K <- length(cur$w); lw <- function(z) log(pmax(z, 1e-300))
  C <- outer(seq_len(K), seq_len(K), Vectorize(function(i, j)
        (prev$m[i] - cur$m[j])^2 + (log(prev$v[i]) - log(cur$v[j]))^2 +
        0.25 * (lw(prev$w[i]) - lw(cur$w[j]))^2))
  C[!is.finite(C)] <- .Machine$double.xmax
  perm <- as.integer(solve_LSAP(C))
  list(w = cur$w[perm], m = cur$m[perm], v = cur$v[perm])
}

## ---- KL against the exact density -----------------------------------------
kl_eval <- function(nu = NULL, w, m, v, n_eval = N_EVAL) {
  rg <- quad_range(nu, tail = KL_TAIL)
  x  <- seq(rg[1], rg[2], length.out = n_eval)
  dx <- (rg[2] - rg[1]) / (n_eval - 1)
  f  <- dens_target(x, nu); f[!is.finite(f)] <- 0; f <- f / (.csum(f) * dx)
  lg <- log_mix(x, w, m, v)
  ok <- f > 1e-300 & is.finite(lg)
  .csum(f[ok] * (log(f[ok]) - lg[ok])) * dx
}

## ---- fit one nu, with adaptive quadrature refinement ----------------------
fit_at_nu <- function(K, nu, prev, n_quad = N_QUAD, n_quad_max = N_QUAD_MAX) {
  qd <- make_quad(nu, n_quad)
  st <- fit_step(K, qd, prev)
  nq <- length(qd$x)
  while (min(sqrt(st$v)) <= 2 * qd$delta * 1.001 && nq < n_quad_max) {
    nq <- 2L * nq; qd <- make_quad(nu, nq)
    st <- fit_step(K, qd, prev)
  }
  st$n_quad <- length(qd$x)
  st
}

## ---- seeding ladder at the first nu: K = 2 .. KMAX ------------------------
seed_ladder <- function(nu, Kmax = KMAX) {
  out <- list(); prev <- NULL
  for (K in 2:Kmax) {
    t0 <- proc.time()[["elapsed"]]
    st <- fit_at_nu(K, nu, prev)
    out[[as.character(K)]] <- st
    cat(sprintf("  ladder nu=%.1f K=%2d  F=%.10f  quad=%5d  %-16s  %5.1fs\n",
                nu, K, st$F, st$n_quad, st$stop_reason,
                proc.time()[["elapsed"]] - t0)); flush.console()
    prev <- list(w = st$w, m = st$m, v = st$v)
  }
  out
}

## ---- walk a sequence of nodes at fixed K ----------------------------------
## `seed` is the converged mixture at grid index `i0`; `idx` is the ordered
## sequence of indices to visit from there (ascending or descending).  Returns
## the fits at `idx` only -- the seed node is recorded by the caller.
walk_nodes <- function(K, seed, i0, idx, nu_grid = NU_GRID) {
  n <- length(idx)
  if (n == 0L) return(list(idx = integer(0)))
  W <- matrix(NA_real_, n, K); M <- matrix(NA_real_, n, K); V <- matrix(NA_real_, n, K)
  KLv <- numeric(n); NQ <- integer(n); SR <- character(n)
  prev <- list(w = seed$w, m = seed$m, v = seed$v); nu_prev <- nu_grid[i0]
  t0 <- proc.time()[["elapsed"]]
  for (j in seq_len(n)) {
    i <- idx[j]
    ## SCALE PREDICTOR.  Part of the change from one node to the next is a pure
    ## rescaling that is known in closed form (sd_target above), so apply it to
    ## the warm start before the solver runs: the optimiser then only has to fix
    ## the residual SHAPE change instead of rediscovering the rescaling.  This
    ## is the predictor half of a predictor-corrector continuation.  It matters
    ## at small nu, where sd ~ 2/sqrt(nu) (r = 0.9875 per node at nu = 0.1), and
    ## is inert at large nu, where sd -> 2 (r = 0.9999 at nu = 20).
    r    <- sd_target(nu_grid[i]) / sd_target(nu_prev)
    warm <- list(w = prev$w, m = prev$m * r, v = prev$v * r^2)
    st <- fit_at_nu(K, nu_grid[i], warm)
    ## match against the PREDICTOR, not the raw previous node: it is the point
    ## the solver actually started from, so it is the closer reference.
    tr <- track(warm, st)
    W[j, ] <- tr$w; M[j, ] <- tr$m; V[j, ] <- tr$v
    KLv[j] <- kl_eval(nu_grid[i], tr$w, tr$m, tr$v)
    NQ[j] <- st$n_quad; SR[j] <- st$stop_reason
    prev <- list(w = tr$w, m = tr$m, v = tr$v); nu_prev <- nu_grid[i]
    if (j %% 25 == 0 || j == n)
      cat(sprintf("  K=%2d  nu=%8.4f  (%3d/%3d)  KL=%.4e  quad=%5d  %6.1fs\n",
                  K, nu_grid[i], j, n, KLv[j], NQ[j],
                  proc.time()[["elapsed"]] - t0)); flush.console()
  }
  list(idx = idx, w = W, m = M, v = V, KL = KLv, n_quad = NQ, stop_reason = SR)
}

## ---- assemble the table: the starting node plus the walk ------------------
.assemble <- function(K, i0, start, walk, nu_grid = NU_GRID) {
  n <- length(nu_grid)
  W <- matrix(NA_real_, n, K); M <- matrix(NA_real_, n, K); V <- matrix(NA_real_, n, K)
  KLv <- rep(NA_real_, n); NQ <- rep(NA_integer_, n); SR <- rep(NA_character_, n)
  W[i0, ] <- start$w; M[i0, ] <- start$m; V[i0, ] <- start$v
  KLv[i0] <- kl_eval(nu_grid[i0], start$w, start$m, start$v)
  NQ[i0] <- start$n_quad; SR[i0] <- start$stop_reason
  if (length(walk$idx)) {
    W[walk$idx, ] <- walk$w; M[walk$idx, ] <- walk$m; V[walk$idx, ] <- walk$v
    KLv[walk$idx] <- walk$KL; NQ[walk$idx] <- walk$n_quad
    SR[walk$idx] <- walk$stop_reason
  }
  list(K = K, nu = nu_grid, w = W, m = M, v = V, KL = KLv,
       nats = TOBS * KLv, n_quad = NQ, stop_reason = SR)
}

## ---- run the walks, sequentially or in parallel ---------------------------
## The walks for different K are independent of one another; only the nodes
## WITHIN a walk are sequential (that is what continuation means), so a single
## walk cannot be split and the speed-up saturates at one worker per K.
## Nothing here draws random numbers, so forked workers need no seed handling.
.run_jobs <- function(jobs, n_workers = N_WORKERS) {
  if (n_workers <= 1L || length(jobs) <= 1L) return(lapply(jobs, function(f) f()))
  parallel::mclapply(jobs, function(f) f(), mc.cores = min(n_workers, length(jobs)),
                     mc.preschedule = FALSE)
}
