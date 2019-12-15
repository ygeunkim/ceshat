#' Inverting 1 - conditional CDF
#'
#' @description
#' This function computes true CVaR given true conditional CDF.
#' @param p Upper tail probability. \code{0.95}, by default.
#' @param cdf conditional cdf given x. Function of \code{(y, x)}
#' @param x conditional information variable \code{x}
#' @param lower lower value to be primarily searched
#' @param upper upper value to be primarily searched
#' @return
#' CVaR value for each \code{x}
#' @details
#' Note that
#' \deqn{\hat{nu}_p(x) = \hat{S}_c^{-1}(p \mid x)}
#' where
#' \deqn{\hat{S}(y \mid x)_c(y \mid x) = 1 - \hat{F}_c(y \mid x)}
#' Inverting can be done by finding
#' \deqn{\inf \{ y : F(y \mid x) \ge 1 - p \}}
#' When finding the minimum, firstly find in \code{(lower, upper)}.
#' But the answer can be out of this bound.
#' These arguments just set the primary interval, not the final one.
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @export
invert_cvar <- function(p = .95, cdf, x, lower = -1, upper = 1) {
  Vectorize(invert_cdf, vectorize.args = "x")(p, cdf, x, lower, upper)
}

invert_cdf <- function(p = .95, cdf, x, lower = -1, upper = 1) {
  y_grid <- seq(lower, upper, by = .01)
  cand <- explore_grid(y_grid, p, cdf, x)
  if (length(cand) > 0) {
    if (min(cand) > lower) {
      return(min(cand))
    } else {
      y_grid <- seq(lower - 5, upper, by = .01)
      return(min(explore_grid(y_grid, p, cdf, x)))
    }
  } else {
    y_grid <- seq(lower, upper + 5, by = .01)
    return(min(explore_grid(y_grid, p, cdf, x)))
  }
}

explore_grid <- function(grids, p, cdf, x) {
  loss <- cdf(grids, x)
  grids[loss >= 1 - p]
  # min(grids[loss >= 1 - p])
}

#' Plug-in Method
#'
#' @description
#' This function computes true CES given true CVaR and conditional pdf.
#' @param p Upper tail probability. \code{0.95}, by default.
#' @param pdf conditional pdf given x. Function of \code{(y, x)}
#' @param cdf conditional cdf or cvar
#' @param x conditional information variable \code{x}
#' @param lower lower value to be searched
#' @param upper upper value to be searched
#' @return
#' CES function of \code{x}
#' @details
#' Note that
#' \deqn{\hat{\mu}_p(x) = \frac{1}{p} \sum_{t = 1}^n W_{c,t}(x, h) \left[ Y_t \bar{G}_{h_0} (\hat{\nu}_p (x) - Y_t) + h_0 G_{1, h_0} (\hat{\nu}_p (x) - Y_t) \right]}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @importFrom stats integrate
#' @export
plugin_ces <- function(p = .95, pdf, cdf, x, lower = -2, upper = 2) {
  Vectorize(plugin_pdf, vectorize.args = "x")(p, pdf, cdf, x, lower, upper)
}

plugin_pdf <- function(p = .95, pdf, cdf, x, lower = -2, upper = 2) {
  y_grid <- seq(lower, upper, by = .01)
  CVAR <- invert_cdf(p, cdf, x)
  integrate(
    function(y) {
      y * pdf(y, x)
    },
    lower = CVAR,
    upper = Inf
  )$value / p
}
