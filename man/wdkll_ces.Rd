% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cai_ces.R
\name{wdkll_ces}
\alias{wdkll_ces}
\title{Weighted Double Kernel Local Linear Estimation of Conditional Expected Shortfall}
\usage{
wdkll_ces(
  formula,
  data,
  prob = 0.95,
  nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"),
  nw_h,
  pdf_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"),
  h0,
  init = 0,
  eps = 1e-05,
  iter = 1000,
  lower_invert = -3,
  upper_invert = 3
)
}
\arguments{
\item{formula}{an object class \link[stats]{formula}.}

\item{data}{an optional data to be used.}

\item{prob}{upper tail probability for VaR}

\item{nw_kernel}{Kernel for weighted nadaraya watson}

\item{nw_h}{Bandwidth for WNW. If not specified, use the asymptotic optimal.}

\item{pdf_kernel}{Kernel for initial estimate of conditinal pdf}

\item{h0}{Bandwidth for pdf kernel. If not specified, use 0.1 times of asymptotic optimal for \code{nw_h}.}

\item{init}{initial value for finding lambda}

\item{eps}{small value}

\item{iter}{maximum iteration when finding lambda}

\item{lower_invert}{lower y when inverting the cdf}

\item{upper_invert}{upper y when inverting the cdf}
}
\description{
WDKLL estimator of CES
}
\details{
Plugging-in in methods gives
\deqn{\hat{\mu}_p(x) = \frac{1}{p} \sum_{t = 1}^n W_{c,t}(x, h) \left[ Y_t \bar{G}_{h_0} (\hat{\nu}_p (x) - Y_t) + h_0 G_{1, h_0} (\hat{\nu}_p (x) - Y_t) \right]}
}
\references{
Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
}
