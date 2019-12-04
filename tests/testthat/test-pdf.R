context("Conditional density")

test_that(
  "WDKLL estimation of Conditional PDF",
  {
    X <- data.frame(x = 1:10, y = rnorm(10))
    newX <- data.frame(x = 11:12, y = rnorm(2))
    fhat <-
      wdkll_pdf(
        y ~ x,
        data = X,
        newdata = newX,
        nw_kernel = kernel_gauss, nw_h = 1,
        pdf_kernel = kernel_gauss, h0 = 1
      )
    expect_type(fhat, "double")
    expect_length(fhat, 2)
  }
)
