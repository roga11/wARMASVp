// ============================================================================ //
// wARMASVp: Core Estimation Functions
// Authors: Nazmul Ahsan, Jean-Marie Dufour, Gabriel Rodriguez Rondon
// ============================================================================ //
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declarations from utils_cpp.cpp
arma::mat acov_g(arma::mat y, int k);
double kendall_corr(arma::vec x, arma::vec y);

// [[Rcpp::export]]
List svpCpp(arma::vec u, int p, int J, bool trunc_lev, double del, int rho_type = 1, bool wDecay = false) {
  Rcpp::Environment gsignal("package:gsignal");
  Rcpp::Function polyR = gsignal["poly"];

  // ----- get y*
  int N = u.n_elem;
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
  // Use decaying weights
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

  // ----- compute rho (delta in paper)
  arma::vec yabs = abs(u);
  double muu = sum(u.subvec(0, N - 2)) / (N - 1);
  double mua = sum(yabs.subvec(1, N - 1)) / (N - 1);
  double EH = 0.0;
  if (rho_type == 2) {
    double EH_kentmp = kendall_corr(yabs.subvec(1, N - 1) - mua, u.subvec(0, N - 2) - muu);
    arma::vec VYA = (1.0 / (N - 2)) * trans(yabs.subvec(1, N - 1) - mua) * (yabs.subvec(1, N - 1) - mua);
    arma::vec VYN = (1.0 / (N - 2)) * trans(u.subvec(0, N - 2) - muu) * (u.subvec(0, N - 2) - muu);
    EH = EH_kentmp * sqrt(VYA(0) * VYN(0));
  } else if (rho_type == 1) {
    EH = ((1.0 / (N - 2)) * as_scalar(trans(yabs.subvec(1, N - 1) - mua) * (u.subvec(0, N - 2) - muu)));
  }
  // gamma tilde
  double gammatilde;
  if (p == 1) {
    gammatilde = pow(s_vh_ols, 2) / (1 - phiB(0));
  } else if (p == 2) {
    gammatilde = pow(s_vh_ols, 2) / ((1 - phiB(0) - phiB(1)) * (1 + phiB(1)));
  } else if (p == 3) {
    gammatilde = ((1 - phiB(2)) * pow(s_vh_ols, 2)) / ((1 - phiB(0) - phiB(1) - phiB(2)) * (1 + phiB(0) * phiB(2) + phiB(1) - pow(phiB(2), 2)));
  } else if (p == 4) {
    gammatilde = ((1 - phiB(2) - phiB(0) * phiB(3) - pow(phiB(3), 2)) * pow(s_vh_ols, 2)) / ((1 - phiB(0) - phiB(1) - phiB(2) - phiB(3)) * (1 + phiB(0) * (phiB(2) + phiB(1)) + phiB(3) + 2 * phiB(1) * phiB(3) + pow(phiB(0), 2) * phiB(3) + phiB(1) * pow(phiB(3), 2) - pow(phiB(2), 2) - pow(phiB(3), 2) - pow(phiB(3), 3) - phiB(0) * phiB(2) * phiB(3)));
  } else {
    Rcpp::warning("Order 'p' must be <=4 for gammatilde computation.");
    gammatilde = NA_REAL;
  }
  // compute rho (delta)
  double rho = (sqrt(2 * arma::datum::pi) * EH) / (s_vh_ols * pow(s_yh, 2)) * exp(-0.25 * gammatilde);
  if (trunc_lev) {
    if (rho < -1) {
      rho = -1;
    } else if (rho > 1) {
      rho = 1;
    }
  }
  // ----- Organize output
  List para = List::create(
    Named("mu") = mu,
    Named("phi") = phiB,
    Named("sigv") = s_vh_ols,
    Named("sigy") = s_yh,
    Named("rho") = rho,
    Named("gammatilde") = gammatilde,
    Named("nonstationary_ind") = (indices.n_elem > 0)
  );
  return para;
}


// [[Rcpp::export]]
List svpCpp_nolev(arma::vec u, int p, int J, double del, bool wDecay = false) {
  Rcpp::Environment gsignal("package:gsignal");
  Rcpp::Function polyR = gsignal["poly"];

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


// [[Rcpp::export]]
List svpCpp_h_ahead(arma::vec u, int p, int J, int h, double del, bool wDecay = false) {
  Rcpp::Environment gsignal("package:gsignal");
  Rcpp::Function polyR = gsignal["poly"];

  int K_min = std::max(p, h + 1);

  // ----- get y*
  arma::vec u2 = square(u) + del;
  double mu = mean(log(u2));
  arma::mat y = log(u2) - mu;
  int T = y.n_elem;

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
        int lag = jj + K_min + i - j + 1;
        xMv(j) = acov_g(y, lag)(0, 0);
      }
      xM.row(i) = xMv;
      int yLag = jj + K_min + i - h + 1;
      yV(i) = acov_g(y, yLag)(0, 0);
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

  // ----- compute sigma_v
  double sigma2_eps = arma::datum::pi * arma::datum::pi / 2.0;
  double phi_sq_sum = arma::dot(phiB, phiB);

  arma::vec resids(T - h, arma::fill::zeros);
  for (int t = h; t < T; ++t) {
    double pred = 0.0;
    for (int j = 0; j < p; ++j) {
      if ((t - j - 1) >= 0) {
        pred += phiB(j) * y(t - j - 1);
      }
    }
    resids(t - h) = y(t) - pred;
  }
  double Vhat = arma::mean(arma::square(resids));
  double s_vh2_ols = Vhat - sigma2_eps * (1.0 + phi_sq_sum);
  double s_vh_ols = (s_vh2_ols < 0) ? sqrt(sqrt(pow(s_vh2_ols, 2))) : sqrt(s_vh2_ols);

  List para = List::create(
    Named("mu") = mu,
    Named("phi") = phiB,
    Named("sigv") = s_vh_ols,
    Named("sigy") = s_yh
  );
  return para;
}
