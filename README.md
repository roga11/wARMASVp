# wARMASVp

Winsorized ARMA Estimation for Higher-Order Stochastic Volatility Models

## Overview

**wARMASVp** provides estimation, simulation, hypothesis testing, and forecasting for univariate higher-order stochastic volatility SV(p) models. It supports Gaussian, Student-t, and Generalized Error Distribution (GED) innovations, with optional leverage effects.

The estimation method is based on closed-form Winsorized ARMA-SV (W-ARMA-SV) moment-based estimators that avoid numerical optimization, making them fast and reliable.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("roga11/wARMASVp")
```

## Features

- **Estimation**: SV(p) models with Gaussian, Student-t, or GED errors via `svp()`
- **Leverage effects**: Asymmetric volatility estimation for Gaussian SV(p)
- **Simulation**: Generate SV(p) data with `sim_svp()`
- **Hypothesis testing**: LMC and MMC procedures for autoregressive order, leverage, and heavy tails
- **Forecasting**: Kalman filter-based h-step-ahead volatility forecasts via `forecast_svp()`
- **Standard errors**: Simulation-based confidence intervals via `svpSE()`

## Quick Start

```r
library(wARMASVp)

# Simulate Gaussian SV(1)
y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)

# Estimate
fit <- svp(y, p = 1)
summary(fit)

# Standard errors
se <- svpSE(fit, n_sim = 99)
se$CI

# Forecast
fc <- forecast_svp(fit, H = 10)
plot(fc)
```

## References

- Ahsan, M. N. and Dufour, J.-M. (2021). Simple estimators and inference for higher-order stochastic volatility models. *Journal of Econometrics*, 224(1), 181-197. [doi:10.1016/j.jeconom.2021.03.008](https://doi.org/10.1016/j.jeconom.2021.03.008)

- Ahsan, M. N., Dufour, J.-M., and Rodriguez-Rondon, G. (2025). Estimation and inference for higher-order stochastic volatility models with leverage. *Journal of Time Series Analysis*, 46(6), 1064-1084. [doi:10.1111/jtsa.12851](https://doi.org/10.1111/jtsa.12851)

- Ahsan, M. N., Dufour, J.-M., and Rodriguez-Rondon, G. (2026). Estimation and inference for stochastic volatility models with heavy-tailed distributions. Bank of Canada Staff Working Paper 2026-8. [doi:10.34989/swp-2026-8](https://doi.org/10.34989/swp-2026-8)

## License

GPL (>= 3)
