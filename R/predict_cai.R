#' Predict method for WDKLL cvar
#'
#' @description
#' WDKLL values for CVar
#' @param object Object of class from \code{\link{wdkll_cvar}}
#' @param newx x to predict. Unless specified, use the \code{data}.
#' @return
#' CVaR given \code{x}
#' @details
#' CVaR can be earned by inverting the CDF.
#' \deqn{\hat{nu}_p(x) = \hat{S}_c^{-1}(p \mid x)}
#' where
#' \deqn{\hat{S}(y \mid x)_c(y \mid x) = 1 - \hat{F}_c(y \mid x)}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @export
predict.cvar <- function(object, newx) {
  if (missing(newx)) newx <- object$xt
  cvar <- Vectorize(object$cvar, vectorize.args = "x")
  cvar(newx)
}

#' Predict method for WDKLL ces
#'
#' @description
#' WDKLL values for CES
#' @param object Object of class from \code{\link{wdkll_ces}}
#' @param newx x to predict. Unless specified, use the \code{data}.
#' @details
#' Plugging-in in methods gives
#' \deqn{\hat{\mu}_p(x) = \frac{1}{p} \sum_{t = 1}^n W_{c,t}(x, h) \left[ Y_t \bar{G}_{h_0} (\hat{\nu}_p (x) - Y_t) + h_0 G_{1, h_0} (\hat{\nu}_p (x) - Y_t) \right]}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @importFrom stats integrate uniroot
#' @export
predict.ces <- function(object, newx) {
  sapply(
    newx,
    function(x) {
      predict_ces(object, x)
    }
  )
}

predict_ces <- function(object, newx) {
  cvar_fit <- object$cvar
  xt <- cvar_fit$xt
  yt <- object$yt
  prob <- cvar_fit$right_tail
  nw_kernel <- cvar_fit$kernel[1]
  nw_h <- cvar_fit$bandwidth[1]
  pdf_kernel <- cvar_fit$kernel[2]
  h0 <- cvar_fit$bandwidth[2]
  init <- object$newton_param[1]
  eps <- object$newton_param[2]
  iter <- object$newton_param[3]
  cvar <- predict(cvar_fit, newx)
  gh0 <- compute_gh(cvar - yt, pdf_kernel, h0)
  g1h <- function(x) {
    x * compute_kernel(x, pdf_kernel, h0)
  }
  g1 <-
    sapply(
      yt,
      function(y) {
        integrate(
          g1h,
          lower = cvar - y,
          upper = Inf,
          rel.tol = eps,
          abs.tol = eps
        )$value
      }
    )
  if (missing(newx)) newx <- xt
  pt <- find_weight(xt, newx, nw_kernel, nw_h, init, eps, iter)
  wh <- compute_kernel(newx - xt, nw_kernel, nw_h)
  wct <- pt * wh / sum(pt * wh)
  sum( wct * ( yt * (1 - gh0) + h0 * g1 )) / prob
}
