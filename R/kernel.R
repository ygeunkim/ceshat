#' Symmetric kernel function
#'
#' @description
#' Kernel function with bandwidth
#' @param x point where the kernel function is evaluated
#' @param type which kernel to use? \code{c("Gaussian", "Epanechinikov", "Tricube", "Boxcar")}.
#' @param h bandwidth
#' @return
#' Kernel function value given \code{x}
#' @details
#' Gaussian kernel function.
#' \deqn{K_h(x) = \frac{1}{\sqrt{2\pi h}} \exp( - \frac{x^2}{2h})}
#' @importFrom stats dnorm
#' @export
compute_kernel <- function(x, type = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), h) {
  type <- match.arg(type)
  if (type == "Gaussian") {
    kernel_gauss(x, h)
  } else if (type == "Epanechinikov") {
    kernel_epa(x, h)
  } else if (type == "Tricube") {
    kernel_tricube(x, h)
  } else if (type == "Boxcar") {
    kernel_boxcar(x, h)
  }
}

kernel_gauss <- function(x, h) {
  dnorm(x, sd = h)
}

kernel_epa <- function(x, h) {
  3 * (1 - x / h)^2 / (4 * h) * (x >= 0)
}

kernel_tricube <- function(x, h) {
  70 * (1 - abs(x / h)^3)^3 / (81 * h) * (x >= 0)
}

kernel_boxcar <- function(x, h) {
  (x >= 0) / (2 * h)
}

#' Distribution function of kernel
#'
#' @description
#' Distribution fuction of \code{\link{compute_kernel}}
#' @param x point where the kernel function is evaluated
#' @param type which kernel to use? \code{c("Gaussian", "Epanechinikov", "Tricube", "Boxcar")}.
#' @param h bandwidth
#' @return
#' Distribution function value
#' @importFrom stats pnorm
#' @export
compute_gh <- function(x, type = c("Gaussian", "Epanechinikov", "Tricube", "Boxcar"), h) {
  type <- match.arg(type)
  x <- x / h
  if (type == "Gaussian") {
    # pnorm(x, sd = h)
    pnorm(x)
  } else if (type == "Epanechinikov") {
    # x * (3 * h^2 - 3 * h * x + x^2) / (4 * h^3) * (x >= 0)
    x <- x / h
    x * (x^2 - 3 * x + 3) / 4 * (x >= 0)
  } else if (type == "Tricube") {
    # x * (140 * h^9 - 150 * h^6 * x^3 + 60 * h^3 * x^6 - x^9) / (162 * h^10) * (x >= 0)
    x * (-14 * x^9 + 60 * x^6 - 105 * x^3 + 140) / 162 * (x >= 0)
  } else if (type == "Boxcar") {
    x / 2 * (x >= 0)
  }
}

