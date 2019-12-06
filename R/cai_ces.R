#' Weighted Double Kernel Local Linear Estimation of Coditional CDF
#'
#' @description
#' This function estimates conditional CDF
#' using WDKLL method.
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
#' @param pdf_kernel Kernel for initial estimate of conditinal pdf
#' @param h0 Bandwidth for pdf kernel
#' @param init initial value for finding lambda
#' @param eps small value
#' @param iter maximum iteration when finding lambda
#' @return
#' Conditional CDF function with argument \code{y} and \code{x}
#' @examples
#' X <- data.frame(x = 1:10, y = rnorm(10))
#' wdkll_cdf(y ~ x, X, nw_h = 1, h0 = 1)(1, 1)
#' wdkll_cdf(y ~ x, X, nw_h = 1, h0 = 1)(2, 1)
#' @details
#' \deqn{\hat{F}_c(y \mid x) = \sum_{t = 1}^n W_{c,t} G_{h_0} (y - Y_t)}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @import dplyr
#' @export
wdkll_cdf <- function(formula, data,
                      nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                      pdf_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), h0,
                      init = 0, eps = 1e-5, iter = 1000) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  nw_kernel <- match.arg(nw_kernel)
  pdf_kernel <- match.arg(pdf_kernel)
  Vectorize(
    function(y, x) {
      wdkll_fit_cdf(xt, yt, x, y, nw_kernel, nw_h, pdf_kernel, h0, init, eps, iter)
    },
    vectorize.args = "y"
  )
}

wdkll_fit_cdf <- function(xt, yt, x, y, nw_kernel, nw_h, pdf_kernel, h0, init = 0, eps = 1e-5, iter = 1000) {
  # distribution function of K
  gh0 <- compute_gh(y - yt, pdf_kernel, h0)
  # WNW
  pt <- find_weight(xt, x, nw_kernel, nw_h, init, eps, iter)
  wh <- compute_kernel(x - xt, nw_kernel, nw_h)
  wct <- pt * wh / sum(pt * wh)
  # WDKLL
  sum( wct * gh0 )
}

#' Weighted Double Kernel Local Linear Estimation of Conditional Value at Risk
#'
#' @description
#' WDKLL estimator of CVaR
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param prob upper tail probability for VaR
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
#' @param pdf_kernel Kernel for initial estimate of conditinal pdf
#' @param h0 Bandwidth for pdf kernel
#' @param init initial value for finding lambda
#' @param eps small value
#' @param iter maximum iteration when finding lambda
#' @param lower_invert lower y when inverting the cdf
#' @param upper_invert upper y when inverting the cdf
#' @return
#' CVaR given \code{x}
#' @details
#' CVaR can be earned by inverting the CDF.
#' \deqn{\hat{nu}_p(x) = \hat{S}_c^{-1}(p \mid x)}
#' where
#' \deqn{\hat{S}(y \mid x)_c(y \mid x) = 1 - \hat{F}_c(y \mid x)}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @import dplyr
#' @import tibble
#' @importFrom tidyr expand_grid
#' @export
wdkll_cvar <- function(formula, data, prob = .95,
                       nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                       pdf_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), h0,
                       init = 0, eps = 1e-5, iter = 1000,
                       lower_invert = -50, upper_invert = 50) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  cvar <- function(x) {
    uniroot(
      function(y) {
        1 -
          wdkll_cdf(formula, data, nw_kernel, nw_h, pdf_kernel, h0, init, eps, iter)(y, x) -
          prob
      },
      lower = lower_invert,
      upper = upper_invert
    )$root
  }
  result <- list(cvar = cvar)
  result$right_tail <- prob
  result$kerel <- c(nw_kernel, pdf_kernel)
  result$bandwidth <- c(nw_h, h0)
  result$xt <- xt
  class(result) <- "cvar"
  result
}

#' Weighted Double Kernel Local Linear Estimation of Conditional Expected Shortfall
#'
#' @description
#' WDKLL estimator of CES
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param prob upper tail probability for VaR
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
#' @param pdf_kernel Kernel for initial estimate of conditinal pdf
#' @param h0 Bandwidth for pdf kernel
#' @param init initial value for finding lambda
#' @param eps small value
#' @param iter maximum iteration when finding lambda
#' @param lower_invert lower y when inverting the cdf
#' @param upper_invert upper y when inverting the cdf
#' @details
#' Plugging-in in methods gives
#' \deqn{\hat{\mu}_p(x) = \frac{1}{p} \sum_{t = 1}^n W_{c,t}(x, h) \left[ Y_t \bar{G}_{h_0} (\hat{\nu}_p (x) - Y_t) + h_0 G_{1, h_0} (\hat{\nu}_p (x) - Y_t) \right]}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @import dplyr
#' @importFrom stats integrate uniroot
#' @export
wdkll_ces <- function(formula, data, prob = .95,
                      nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                      pdf_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), h0,
                      init = 0, eps = 1e-5, iter = 1000,
                      lower_invert = -50, upper_invert = 50) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  cvar_fit <- wdkll_cvar(formula, data, prob, nw_kernel, nw_h, pdf_kernel, h0, init, eps, iter, lower_invert, upper_invert)
  result <- list(cvar = cvar_fit)
  result$yt <- yt
  result$newton_param <- c(init, eps, iter)
  class(result) <- "ces"
  result
}

# CVaR class-------------------------------

#' `cvar` class
#'
#' @description
#' The \code{cvar} class is a result of \code{\link{wdkll_cvar}}.
#' @name cvar-class
#' @rdname cvar-class
#' @aliases cvar cvar-class
#' @importFrom methods setOldClass
#' @exportClass cvar
setOldClass("cvar")

# CES class--------------------------------

#' `ces` class
#'
#' @description
#' The \code{ces} class is a result of \code{\link{wdkll_ces}}.
#' @name ces-class
#' @rdname ces-class
#' @aliases ces ces-class
#' @importFrom methods setOldClass
#' @exportClass ces
setOldClass("ces")
