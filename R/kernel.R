#' Gaussian kernel function
#'
#' @description
#' Gaussian kernel function with bandwidth
#' @param x point where the kernel function is evaluated
#' @param h bandwidth
#' @return
#' Gaussian function value given \code{x}
#' @details
#' \code{h} is the standard deviation.
#' \deqn{K_h(x) = \frac{1}{\sqrt{2\pi h}} \exp( - \frac{x^2}{2h})}
#' @importFrom stats dnorm
#' @export
kernel_gauss <- function(x, h) {
  dnorm(x, sd = h)
}
