// ============================================================================ //
// wARMASVp: Root-Finding for Distribution Parameter Estimation
// Authors: Nazmul Ahsan, Jean-Marie Dufour, Gabriel Rodriguez Rondon
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declaration from utils_cpp.cpp
double qged_std_cpp(double p, double nu);

// --------------------------------------------------------------------------- //
// Brent's method for root-finding (Netlib public domain implementation)
// Finds x in [a,b] such that f(x) = 0, where f(a) and f(b) have opposite signs.
// tol: convergence tolerance on x
// maxiter: maximum iterations
// --------------------------------------------------------------------------- //
// This is a templated Brent solver that takes a C++ callable (lambda/functor)
template<typename F>
static double brent_zero(F f, double a, double b, double tol, int maxiter) {
  double fa = f(a);
  double fb = f(b);
  // If same sign, return the endpoint closer to zero
  if (fa * fb > 0) {
    return (std::fabs(fa) < std::fabs(fb)) ? a : b;
  }

  double c = a, fc = fa;
  double d = b - a, e = d;

  for (int iter = 0; iter < maxiter; iter++) {
    if (fb * fc > 0) {
      c = a; fc = fa;
      d = b - a; e = d;
    }
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }

    double tol1 = 2.0 * std::numeric_limits<double>::epsilon() * std::fabs(b) + 0.5 * tol;
    double m = 0.5 * (c - b);

    if (std::fabs(m) <= tol1 || fb == 0.0) {
      return b;
    }

    if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      // Try inverse quadratic interpolation
      double s = fb / fa;
      double p_val, q;
      if (a == c) {
        p_val = 2.0 * m * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        double r = fb / fc;
        p_val = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p_val > 0) q = -q;
      else p_val = -p_val;

      if (2.0 * p_val < std::min(3.0 * m * q - std::fabs(tol1 * q), std::fabs(e * q))) {
        e = d;
        d = p_val / q;
      } else {
        d = m; e = m;
      }
    } else {
      d = m; e = m;
    }

    a = b; fa = fb;
    if (std::fabs(d) > tol1) {
      b += d;
    } else {
      b += (m > 0 ? tol1 : -tol1);
    }
    fb = f(b);
  }
  return b;
}


// =========================================================================== //
// Student-t: solve psigamma(nu/2, 1) = se2b for nu
// =========================================================================== //

// [[Rcpp::export]]
double find_nu_t_cpp(double se2b, double nu_lower = 2.01, double nu_upper = 500.0,
                      bool logNu = false, double tol = 1e-6, int maxiter = 1000) {
  // Boundary checks (same as R)
  if (se2b <= R::psigamma(nu_upper / 2.0, 1)) {
    return nu_upper;
  }
  if (se2b >= R::psigamma(nu_lower / 2.0, 1)) {
    return nu_lower;
  }

  if (logNu) {
    double log_lower = std::log(nu_lower);
    double log_upper = std::log(nu_upper);
    auto f_log = [se2b](double logx) -> double {
      return R::psigamma(std::exp(logx) / 2.0, 1) - se2b;
    };
    double log_root = brent_zero(f_log, log_lower, log_upper, tol, maxiter);
    return std::exp(log_root);
  } else {
    auto f = [se2b](double x) -> double {
      return R::psigamma(x / 2.0, 1) - se2b;
    };
    return brent_zero(f, nu_lower, nu_upper, tol, maxiter);
  }
}


// =========================================================================== //
// GED: solve (2/nu)^2 * psigamma(1/nu, 1) = se2 for nu
// =========================================================================== //

// [[Rcpp::export]]
double find_nu_ged_cpp(double se2, double lower, double upper,
                        double tol = 1e-6, int maxiter = 1000) {
  auto f = [se2](double x) -> double {
    return std::pow(2.0 / x, 2) * R::psigamma(1.0 / x, 1) - se2;
  };

  double fl = f(lower);
  double fu = f(upper);

  if (fl * fu < 0) {
    return brent_zero(f, lower, upper, tol, maxiter);
  }

  // Fallback to fixed interval
  double fixed_lower = 0.1;
  double fixed_upper = 20.0;
  double ffl = f(fixed_lower);
  double ffu = f(fixed_upper);

  if (ffl * ffu < 0) {
    return brent_zero(f, fixed_lower, fixed_upper, tol, maxiter);
  }

  // Boundary returns
  if (se2 <= std::pow(2.0 / fixed_upper, 2) * R::psigamma(1.0 / fixed_upper, 1)) {
    return fixed_upper;
  }
  return fixed_lower;
}


// =========================================================================== //
// GED leverage: solve E[g(z + delta*sigv/2)] = target for delta
// Uses Gauss-Hermite quadrature with pre-computed nodes/weights
// =========================================================================== //

// [[Rcpp::export]]
double find_delta_ged_cpp(double target, double sigv, double nu,
                           arma::vec gh_nodes, arma::vec gh_weights,
                           double tol = 1e-8) {
  int n_gh = gh_nodes.n_elem;

  auto f_root = [&](double d) -> double {
    double c_shift = d * sigv / 2.0;
    double E_g = 0.0;
    for (int i = 0; i < n_gh; i++) {
      double p_val = R::pnorm5(gh_nodes(i) + c_shift, 0.0, 1.0, 1, 0);
      E_g += gh_weights(i) * qged_std_cpp(p_val, nu);
    }
    return E_g - target;
  };

  // Check bounds (match R behavior)
  double f_lo, f_hi;
  bool lo_ok = true, hi_ok = true;
  try { f_lo = f_root(-0.999); } catch (...) { lo_ok = false; }
  try { f_hi = f_root(0.999); } catch (...) { hi_ok = false; }

  if (!lo_ok || !hi_ok || f_lo * f_hi >= 0) {
    if (!lo_ok || !hi_ok) {
      return NA_REAL;
    }
    // Return boundary
    return (std::fabs(f_lo) < std::fabs(f_hi)) ? -0.999 : 0.999;
  }

  return brent_zero(f_root, -0.999, 0.999, tol, 1000);
}
