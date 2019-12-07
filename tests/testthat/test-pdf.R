context("Conditional density")

test_that(
  "WDKLL estimation of Conditional PDF",
  {
    X <- data.frame(x = 1:10, y = rnorm(10))
    fhat <-
      wdkll_pdf(
        y ~ x,
        data = X,
        nw_kernel = "Gaussian", nw_h = 1,
        pdf_kernel = "Gaussian", h0 = 1
      )(rnorm(2), 2)
    expect_type(fhat, "double")
    expect_length(fhat, 2)
  }
)
