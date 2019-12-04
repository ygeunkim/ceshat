#' Symmetric kernel function
#'
#' @description
#' Kernel function with bandwidth
#' @param x point where the kernel function is evaluated
#' @param type which kernel to use? \code{c("Gaussian", "Epanechinikov", "Tricube", "Boxcar")}.
#' @param h bandwidth
#' @return
#' Gaussian function value given \code{x}
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

