// ============================================================================ //
// wARMASVp: Winsorized ARMA Estimation for Higher-Order SV Models
// Authors: Nazmul Ahsan, Jean-Marie Dufour, Gabriel Rodriguez Rondon
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double kendall_corr(arma::vec x, arma::vec y) {
  int n = x.n_elem;
  if (n != (int)y.n_elem) {
    Rcpp::stop("Input vectors must have the same length.");
  }
  int concordant = 0;
  int discordant = 0;
  int tied_pairs = 0;
  for (int i = 0; i < (n - 1); ++i) {
    for (int j = i + 1; j < n; ++j) {
      double diff_x = x(i) - x(j);
      double diff_y = y(i) - y(j);
      if (diff_x == 0 && diff_y == 0) {
        tied_pairs++;
      } else {
        concordant += (diff_x * diff_y > 0);
        discordant += (diff_x * diff_y < 0);
      }
    }
  }
  double adj_concordant = concordant + 0.5 * tied_pairs;
  double adj_discordant = discordant + 0.5 * tied_pairs;
  if (adj_concordant + adj_discordant == 0) {
    return 0.0;
  }
  double tau = (adj_concordant - adj_discordant) / sqrt((adj_concordant + adj_discordant) * n * (n - 1) / 2);
  return tau;
}

arma::mat acov_g(arma::mat y, int k) {
  int Tsize = y.n_rows;
  arma::mat gamma = (1.0 / (Tsize - k)) * trans(y.rows(k, Tsize - 1)) * y.rows(0, Tsize - k - 1);
  return gamma;
}

// Quantile function for standardized GED(nu) with Var=1
// Uses the relationship: F_GED^{-1}(p) = a * qgamma(2*(p-0.5), 1/nu, 1)^{1/nu}
// where a = sqrt(Gamma(1/nu) / Gamma(3/nu)).
// Input p is clamped to [1e-15, 1-1e-15].
double qged_std_cpp(double p, double nu) {
  // Clamp to avoid Inf from qgamma
  if (p < 1e-15) p = 1e-15;
  if (p > 1.0 - 1e-15) p = 1.0 - 1e-15;
  double a = std::exp(0.5 * (R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu)));
  if (p >= 0.5) {
    return a * std::pow(R::qgamma(2.0 * (p - 0.5), 1.0 / nu, 1.0, 1, 0), 1.0 / nu);
  } else {
    return -a * std::pow(R::qgamma(2.0 * (0.5 - p), 1.0 / nu, 1.0, 1, 0), 1.0 / nu);
  }
}

// CDF of standardized GED(nu) with Var=1.  Companion to qged_std_cpp.
// F_GED(u) for u >= 0:  0.5 + 0.5 * F_Gamma((u/a)^nu; 1/nu, 1)
// For u < 0 use symmetry.   a = sqrt(Gamma(1/nu) / Gamma(3/nu)).
// Clamp tail probabilities to [1e-15, 1-1e-15] to match qged_std_cpp.
double pged_std_cpp(double u, double nu) {
  double a = std::exp(0.5 * (R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu)));
  double p;
  if (u >= 0.0) {
    p = 0.5 + 0.5 * R::pgamma(std::pow(u / a, nu), 1.0 / nu, 1.0, 1, 0);
  } else {
    p = 0.5 - 0.5 * R::pgamma(std::pow(-u / a, nu), 1.0 / nu, 1.0, 1, 0);
  }
  if (p < 1e-15)        p = 1e-15;
  if (p > 1.0 - 1e-15)  p = 1.0 - 1e-15;
  return p;
}

// GED copula proxy for the leverage shift and its derivative, SVHT Section 3.3.2:
//   g(u)  = Phi^{-1}( F_GED(u; nu) )            (SVHT eq. 49),
//   g'(u) = f_GED(u; nu) / phi( g(u) ),
// with f_GED(u) = nu / (2 a Gamma(1/nu)) exp(-(|u|/a)^nu) the standardized GED(nu)
// density (Var = 1), a = sqrt(Gamma(1/nu) / Gamma(3/nu)), phi the N(0,1) density.
// Tail-stable: F_GED(u) rounds to 1 in double precision once |u| is large, so the
// inversion uses the UPPER-tail probability of |u| on the log scale,
//   log q = log(1/2) + log Q_Gamma((|u|/a)^nu; 1/nu),   g(|u|) = Phi^{-1}(1 - q),
// via R's qnorm(log_p, lower_tail = FALSE) and pgamma(lower_tail = FALSE, log_p),
// and g' in logs. No probability clamp (unlike pged_std_cpp): g is odd, strictly
// increasing and finite for every finite u, g(u) ~ sqrt(2) (|u|/a)^{nu/2}.
// g' enters the Jacobian A_t of the prediction covariance (SVHT eq. 48b).
void ged_copula_cpp(double u, double nu, double& g, double& dg) {
  const double a     = std::exp(0.5 * (R::lgammafn(1.0 / nu) - R::lgammafn(3.0 / nu)));
  const double x     = std::pow(std::fabs(u) / a, nu);                      // Gamma(1/nu, 1) variate
  const double log_q = std::log(0.5) + R::pgamma(x, 1.0 / nu, 1.0, 0, 1);  // log upper-tail prob of |u|
  const double g_abs = R::qnorm(log_q, 0.0, 1.0, 0, 1);                     // upper-tail N(0,1) quantile
  g  = (u < 0.0) ? -g_abs : g_abs;
  const double log_f = std::log(nu) - std::log(2.0 * a) - R::lgammafn(1.0 / nu) - x;
  dg = std::exp(log_f - R::dnorm(g, 0.0, 1.0, 1));
}

// [[Rcpp::export]]
arma::vec rged_cpp(int n, double mean = 0.0, double sd = 1.0, double nu = 2.0) {
  // Generalized Error Distribution random deviates
  double lambda = std::sqrt(std::pow(2.0, -2.0 / nu) * R::gammafn(1.0 / nu) / R::gammafn(3.0 / nu));
  arma::vec r(n);
  arma::vec g = Rcpp::rgamma(n, 1.0 / nu, 1.0);
  arma::vec u = Rcpp::runif(n);
  for (int i = 0; i < n; ++i) {
    double sign = (u[i] < 0.5) ? -1.0 : 1.0;
    r[i] = lambda * std::pow(2.0 * g[i], 1.0 / nu) * sign;
  }
  return mean + sd * r;
}
