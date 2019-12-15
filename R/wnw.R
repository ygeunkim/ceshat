#' Weighted Nadaraya Watson Estimator of Coditional PDF
#'
#' @description
#' This function estimates conditional pdf
#' using WDKLL method.
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param wt weights for WNW. Computing in prediction step will help efficiency.
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
#' @param init initial value for finding lambda
#' @param eps small value
#' @param iter maximum iteration when finding lambda
#' @return
#' Conditional pdf function of \code{(y, x)}. \code{y} can be a numeric vector.
#' @details
#' Since standalone LL or WNW does not fully satisfy the conditions of cdf,
#' Cai et al (2008) proposed to use WNW in LL scheme.
#' \deqn{\hat{f}_c(y \mid x) = \sum_{t = 1}^n W_{c,t}(w, h) K_{h_0}(y - Y_t)}
#' @import dplyr
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @export
wnw_pdf <- function(formula, data, wt,
                      nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                      init = 0, eps = 1e-5, iter = 1000) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  nw_kernel <- match.arg(nw_kernel)
  Vectorize(
    function(y, x) {
      wnw_fit_pdf(xt, yt, x, y, wt, nw_kernel, nw_h, init, eps, iter)
    },
    vectorize.args = "y"
  )
}

wnw_fit_pdf <- function(xt, yt, x, y, pt, nw_kernel, nw_h, init = 0, eps = 1e-5, iter = 1000) {
  wh <- compute_kernel(x - xt, nw_kernel, nw_h)
  wct <- pt * wh / sum(pt * wh)
  # fhat
  sum( wct * (yt == y) )
}

wnw_fit_pdf2 <- function(xt, yt, wt,
                         nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                         init = 0, eps = 1e-5, iter = 1000) {
  nw_kernel <- match.arg(nw_kernel)
  Vectorize(
    function(y, x) {
      wnw_fit_pdf(xt, yt, x, y, wt, nw_kernel, nw_h, init, eps, iter)
    },
    vectorize.args = "y"
  )
}

#' Weighted Nadaraya Watson Estimator of Coditional CDF
#'
#' @description
#' This function estimates conditional CDF
#' by WNW.
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param wt weights for WNW. Computing in prediction step will help efficiency.
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
#' @param init initial value for finding lambda
#' @param eps small value
#' @param iter maximum iteration when finding lambda
#' @return
#' Conditional CDF function with argument \code{y} and \code{x}
#' @details
#' \deqn{\hat{F}_c(y \mid x) = \sum_{t = 1}^n W_{c,t} I (Y_t \le y)}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @import dplyr
#' @export
wnw_cdf <- function(formula, data, wt,
                      nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                      init = 0, eps = 1e-5, iter = 1000) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  nw_kernel <- match.arg(nw_kernel)
  Vectorize(
    function(y, x) {
      wnw_fit_cdf(xt, yt, x, y, wt, nw_kernel, nw_h, init, eps, iter)
    },
    vectorize.args = "y"
  )
}

wnw_fit_cdf <- function(xt, yt, x, y, pt, nw_kernel, nw_h, init = 0, eps = 1e-5, iter = 1000) {
  wh <- compute_kernel(x - xt, nw_kernel, nw_h)
  wct <- pt * wh / sum(pt * wh)
  sum( wct * (yt <= y) )
}

wnw_cdf2 <- function(xt, yt, wt,
                       nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                       init = 0, eps = 1e-5, iter = 1000) {
  nw_kernel <- match.arg(nw_kernel)
  Vectorize(
    function(y, x) {
      wnw_fit_cdf(xt, yt, x, y, wt, nw_kernel, nw_h, init, eps, iter)
    },
    vectorize.args = "y"
  )
}

# Conditional Value at Risk---------------------------------------------------

#' Weighted Double Kernel Local Linear Estimation of Conditional Value at Risk
#'
#' @description
#' WDKLL estimator of CVaR
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param prob upper tail probability for VaR
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
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
wnw_cvar <- function(formula, data, prob = .95,
                       nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                       init = 0, eps = 1e-5, iter = 1000,
                       lower_invert = -3, upper_invert = 3) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  nw_kernel <- match.arg(nw_kernel)
  if (missing(nw_h)) nw_h <- length(xt)^(-4 / 5)
  result <- list(cvar = c(lower_invert, upper_invert))
  result$right_tail <- prob
  result$kernel <- nw_kernel
  result$bandwidth <- nw_h
  result$yt <- yt
  result$xt <- xt
  result$newton_param <- c(init, eps, iter)
  class(result) <- "nwcvar"
  result
}

#' Predict method for nwcvar
#'
#' @description
#' WDKLL values for CVar
#' @param object Object of class from \code{\link{wdkll_cvar}}
#' @param newx x to predict. Unless specified, use the \code{data}.
#' @param nw Use NW or WNW. \code{TRUE} if NW.
#' @param ... further arguments passed to or from other methods.
#' @return
#' CVaR given \code{x}
#' @details
#' CVaR can be earned by inverting the CDF.
#' \deqn{\hat{nu}_p(x) = \hat{S}_c^{-1}(p \mid x)}
#' where
#' \deqn{\hat{S}(y \mid x)_c(y \mid x) = 1 - \hat{F}_c(y \mid x)}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @export
predict.nwcvar <- function(object, newx, nw = FALSE, ...) {
  if (missing(newx)) newx <- object$xt
  xt <- object$xt
  yt <- object$yt
  prob <- object$right_tail
  nw_kernel <- object$kernel
  nw_h <- object$bandwidth
  init <- object$newton_param[1]
  eps <- object$newton_param[2]
  iter <- object$newton_param[3]
  if (nw) {
    pt <- rep(1, length(xt))
  } else {
    pt <- find_pt(xt, newx, nw_kernel, nw_h, init, eps, iter)
  }
  sapply(
    1:length(newx),
    function(i) {
      predict_nwcvar(object, newx[i], prob, xt, yt, pt[,i], nw_kernel, nw_h, init, eps, iter)
    }
  )
}

predict_nwcvar <- function(object, newx, prob,
                         xt, yt, pt,
                         nw_kernel, nw_h,
                         init, eps, iter) {
  find_cvar <- seq(object$cvar[1], object$cvar[2], by = .01)
  loss <- wnw_cdf2(xt, yt, pt, nw_kernel, nw_h, init, eps, iter)
  cand <- explore_grid(find_cvar, prob, loss, newx)
  if (length(cand) > 0) {
    if (min(cand) > object$cvar[1]) {
      return(min(cand))
    } else {
      find_cvar <- seq(object$cvar[1] - 5, object$cvar[2], by = .01)
      return(min(explore_grid(find_cvar, prob, loss, newx)))
    }
  } else {
    find_cvar <- seq(object$cvar[1], object$cvar[2] + 5, by = .01)
    return(min(explore_grid(find_cvar, prob, loss, newx)))
  }
}

# Conditional Expected Shortfall---------------------------------------------------

#' Weighted Nadaraya Watson Estimator of Conditional Expected Shortfall
#'
#' @description
#' WNW estimator of CES
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param prob upper tail probability for VaR
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW. If not specified, use the asymptotic optimal.
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
#' @importFrom stats integrate
#' @export
wnw_ces <- function(formula, data, prob = .95,
                      nw_kernel = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), nw_h,
                      init = 0, eps = 1e-5, iter = 1000,
                      lower_invert = -3, upper_invert = 3) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  nw_kernel <- match.arg(nw_kernel)
  if (missing(nw_h)) nw_h <- length(xt)^(-4 / 5)
  cvar_fit <- wnw_cvar(formula, data, prob, nw_kernel, nw_h, init, eps, iter, lower_invert, upper_invert)
  result <- list(cvar = cvar_fit)
  class(result) <- "nwces"
  result
}

#' Predict method for nwces
#'
#' @description
#' WNW values for CES
#' @param object Object of class from \code{\link{wdkll_ces}}
#' @param newx x to predict. Unless specified, use the \code{data}.
#' @param nw Use NW or WNW. \code{TRUE} if NW.
#' @param ... further arguments passed to or from other methods.
#' @details
#' Plugging-in in methods gives
#' \deqn{\hat{\mu}_p(x) = \frac{1}{p} \sum_{t = 1}^n W_{c,t}(x, h) \left[ Y_t \bar{G}_{h_0} (\hat{\nu}_p (x) - Y_t) + h_0 G_{1, h_0} (\hat{\nu}_p (x) - Y_t) \right]}
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @importFrom stats integrate
#' @export
predict.nwces <- function(object, newx, nw = FALSE, ...) {
  cvar_fit <- object$cvar
  # cvar <- predict(cvar_fit, newx)
  xt <- cvar_fit$xt
  yt <- cvar_fit$yt
  prob <- cvar_fit$right_tail
  nw_kernel <- cvar_fit$kernel
  nw_h <- cvar_fit$bandwidth
  init <- cvar_fit$newton_param[1]
  eps <- cvar_fit$newton_param[2]
  iter <- cvar_fit$newton_param[3]
  if (nw) {
    pt <- rep(1, length(xt))
  } else {
    pt <- find_pt(xt, newx, nw_kernel, nw_h, init, eps, iter)
  }
  fhat <- wnw_fit_pdf2(xt, yt, pt, nw_kernel, nw_h, init, eps, iter)
  sapply(
    newx,
    function(x) {
      integrate(
        function(y) {
          y * fhat(y, x)
        },
        lower = predict(cvar_fit, x),
        upper = Inf,
        rel.tol = eps,
        abs.tol = eps
      )$value
    }
  )
}

# CVaR class-------------------------------

#' `nwcvar` class
#'
#' @description
#' The \code{nwcvar} class is a result of \code{\link{wnw_cvar}}.
#' @name nwcvar-class
#' @rdname nwcvar-class
#' @aliases nwcvar nwcvar-class
#' @importFrom methods setOldClass
#' @exportClass cvar
setOldClass("nwcvar")

# CES class--------------------------------

#' `nwces` class
#'
#' @description
#' The \code{nwces} class is a result of \code{\link{wnw_ces}}.
#' @name nwces-class
#' @rdname nwces-class
#' @aliases nwces nwces-class
#' @importFrom methods setOldClass
#' @exportClass ces
setOldClass("nwces")

