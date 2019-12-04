#' Weighted Double Kernel Local Linear Estimation of Coditional PDF
#'
#' @description
#' This function estimates conditional pdf
#' using WDKLL method.
#' @param formula an object class \link[stats]{formula}.
#' @param data an optional data to be used.
#' @param newdata data to predict. If \code{NULL}, use the \code{data}.
#' @param nw_kernel Kernel for weighted nadaraya watson
#' @param nw_h Bandwidth for WNW
#' @param pdf_kernel Kernel for initial estimate of conditinal pdf
#' @param h0 Bandwidth for pdf kernel
#' @param init initial value for finding lambda
#' @param eps small value
#' @param iter maximum iteration when finding lambda
#' @return
#' Conditional pdf values for each x
#' @details
#' Since standalone LL or WNW does not fully satisfy the conditions of cdf,
#' Cai et al (2008) proposed to use WNW in LL scheme.
#' \deqn{\hat{f}_c(y \mid x) = \sum_{t = 1}^n W_{c,t}(w, h) K_{h_0}(y - Y_t)}
#' @import dplyr
#' @references Cai, Z., & Wang, X. (2008). \emph{Nonparametric estimation of conditional VaR and expected shortfall}. Journal of Econometrics, 147(1), 120-130.
#' @export
wdkll_pdf <- function(formula, data, newdata = NULL, nw_kernel, nw_h, pdf_kernel, h0, init = 0, eps = 1e-5, iter = 1000) {
  var_name <- find_name(formula)
  yt <- data %>% select(var_name[1]) %>% pull()
  xt <- data %>% select(var_name[2]) %>% pull()
  if (is.null(newdata)) newdata <- data
  if (!is.numeric(newdata)) {
    newy <- newdata %>% select(var_name[1]) %>% pull()
    newx <- newdata %>% select(var_name[2]) %>% pull()
    sapply(
      1:length(newy),
      function(i) {
        wdkll_fit_pdf(xt, yt, newx[i], newy[i], nw_kernel, nw_h, pdf_kernel, h0, init, eps, iter)
      }
    )
  } else {
    newy <- newdata[var_name[1]]
    newx <- newdata[var_name[2]]
    wdkll_fit_pdf(xt, yt, newx, newy, nw_kernel, nw_h, pdf_kernel, h0, init, eps, iter)
  }
}

emp_log <- function(lambda, xt, x, Kh, h) {
  # log-pdf
  log(1 + lambda * (xt - x) * Kh(x - xt, h))
}

deriv_log <- function(lambda, xt, x, Kh, h) {
  # derivative of each log-pdf
  (xt - x) * Kh(x - xt, h) / (1 + lambda * (xt - x) * Kh(x - xt, h))
}

deriv2_log <- function(lambda, xt, x, Kh, h) {
  # second derivative
  -deriv_log(lambda, xt, x, Kh, h)^2
}

find_weight <- function(xdata, x, Kh, h, init_lambda = 0, eps = 1e-5, max_iter = 1000) {
  if (!is.numeric(xdata)) stop("xdata should be numeric vector")
  lambda <- init_lambda
  grad <- 0
  hess <- 0
  for (i in seq_len(max_iter)) {
    grad <- -sum( deriv_log(init_lambda, xdata, x, Kh, h) )
    hess <- -sum( deriv2_log(init_lambda, xdata, x, Kh, h) )
    hess <- ifelse(hess == 0, hess + eps, hess)
    lambda <- init_lambda - grad / hess
    if (abs(init_lambda - lambda) < eps) break
    init_lambda <- lambda
  }
  # pt = 1 / (n * (1 + lambda * (xt - x) * Kh))
  1 / (length(xdata) * (1 + lambda * (xdata - x) * Kh(x - xdata, h)))
}

find_name <- function(formula) {
  formula %>%
    deparse() %>%
    stringr::str_remove_all(pattern = "[:blank:]") %>%
    stringr::str_split(pattern = "\\~|\\+", simplify = TRUE)
}

wdkll_fit_pdf <- function(xt, yt, x, y, nw_kernel, nw_h, pdf_kernel, h0, init = 0, eps = 1e-5, iter = 1000) {
  ystar <- pdf_kernel(y - yt, h0)
  ystar0 <- pdf_kernel(0, h0)
  # pt using ystar
  pt <- find_weight(xt, x, nw_kernel, nw_h, init, eps, iter)
  # Wct
  wct <- pt * nw_kernel(x - xt, nw_h) / sum(pt * nw_kernel(x - xt, nw_h))
  # fhat
  sum( wct * ystar )
}

