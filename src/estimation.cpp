// ============================================================================ //
// wARMASVp: Core Estimation Functions
// Authors: Nazmul Ahsan, Jean-Marie Dufour, Gabriel Rodriguez Rondon
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declaration from utils_cpp.cpp
arma::mat acov_g(arma::mat y, int k);

// [[Rcpp::export]]
List svpCpp_nolev(arma::vec u, int p, int J, double del, bool wDecay = false) {
  Rcpp::Function polyR = Rcpp::Environment::namespace_env("gsignal")["poly"];

  // ----- get y*
  arma::vec u2 = square(u) + del;
  double mu = mean(log(u2));
  arma::mat y = log(u2) - mu;

  // ----- compute sigma_y
  double s_yh = sqrt(exp(mu + 1.2704));

  // ----- compute phi
  arma::mat xMM1(J * p, p, arma::fill::zeros);
  arma::mat yVM1(J * p, 1, arma::fill::zeros);
  arma::mat wDec(J * p, 1, arma::fill::zeros);
  for (int jj = 0; jj < J; ++jj) {
    arma::mat xM(p, p, arma::fill::zeros);
    arma::mat yV(p, 1, arma::fill::zeros);
    for (int i = 0; i < p; ++i) {
      arma::rowvec xMv(p, arma::fill::zeros);
      for (int j = 0; j < p; ++j) {
        int ki = i + p - 1 + jj;
        int k = ki - j + 1;
        xMv(j) = acov_g(y, k)(0, 0);
      }
      xM.row(i) = xMv;
      yV(i) = acov_g(y, p + i + jj + 1)(0, 0);
    }
    xMM1.rows(jj * p, (jj + 1) * p - 1) = xM;
    yVM1.rows(jj * p, (jj + 1) * p - 1) = yV;
    wDec.rows(jj * p, (jj + 1) * p - 1).fill((2.0 / J) * (1.0 - (double(jj + 1) / (J + 1))));
  }
  arma::mat phiB = pinv(trans(xMM1) * xMM1) * trans(xMM1) * yVM1;
  if (wDecay == true) {
    arma::mat rowones(1, p, arma::fill::ones);
    arma::mat wolsDec = arma::sqrt(wDec);
    arma::mat wolsDecX = arma::sqrt(wDec) * rowones;
    phiB = pinv(trans(wolsDecX % xMM1) * (wolsDecX % xMM1)) * trans(wolsDecX % xMM1) * (wolsDec % yVM1);
  }
  // Check stationarity
  double Del = 1 - 0.0001;
  arma::vec pl = arma::join_cols(arma::vec{1}, -phiB);
  arma::cx_vec r = roots(pl);
  arma::uvec indices = find(abs(r) >= 1);
  r.elem(indices) = sign(r.elem(indices)) * Del;
  arma::vec polyc = real(as<arma::vec>(polyR(r)));
  phiB = -polyc.subvec(1, polyc.n_elem - 1);

  // ----- compute sigma_v
  double s_vh2_ols;
  double gamma0 = as_scalar(arma::var(y)) - pow(arma::datum::pi, 2) / 2;
  arma::vec phigB = trans(phiB) * xMM1.submat(0, xMM1.n_cols - 1, p - 1, xMM1.n_cols - 1);
  if (p == 1) {
    s_vh2_ols = (1 - pow(phiB(0), 2)) * gamma0;
  } else {
    s_vh2_ols = gamma0 - phigB(0);
  }
  double s_vh_ols;
  if (s_vh2_ols < 0) {
    s_vh_ols = sqrt(sqrt(pow(s_vh2_ols, 2)));
  } else {
    s_vh_ols = sqrt(s_vh2_ols);
  }

  // ----- Organize output
  List para = List::create(
    Named("mu") = mu,
    Named("phi") = phiB,
    Named("sigv") = s_vh_ols,
    Named("sigy") = s_yh,
    Named("nonstationary_ind") = (indices.n_elem > 0)
  );
  return para;
}


