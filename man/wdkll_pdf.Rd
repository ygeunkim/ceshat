% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cai.R
\name{wdkll_pdf}
\alias{wdkll_pdf}
\title{Weighted Double Kernel Local Linear Estimation of Coditional PDF}
\usage{
wdkll_pdf(
  formula,
  data,
  nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"),
  nw_h,
  pdf_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"),
  h0,
  init = 0,
  eps = 1e-05,
  iter = 1000
)
}
\arguments{
\item{formula}{an object class \link[stats]{formula}.}

\item{data}{an optional data to be used.}

\item{nw_kernel}{Kernel for weighted nadaraya watson}

\item{nw_h}{Bandwidth for WNW}

\item{pdf_kernel}{Kernel for initial estimate of conditinal pdf}

\item{h0}{Bandwidth for pdf kernel}

\item{init}{initial value for finding lambda}

\item{eps}{small value}

\item{iter}{maximum iteration when finding lambda}
}
\value{
Conditional pdf function of \code{(y, x)}. \code{y} can be a numeric vector.
}
\description{
This function estimates conditional pdf
using WDKLL method.
}
\details{
Since standalone LL or WNW does not fully satisfy the conditions of cdf,
Cai et al (2008) proposed to use WNW in LL scheme.
\deqn{\hat{f}_c(y \mid x) = \sum_{t = 1}^n W_{c,t}(w, h) K_{h_0}(y - Y_t)}
}
\references{
Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
}
