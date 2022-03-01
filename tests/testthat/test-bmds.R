library(maotai)

test_that("Bayesian MDS works", {
  data(iris)
  dat = as.matrix(iris[,1:4])
  
  iris2d = maotai::bmds(dat, ndim=2, mc.iter=5, par.step=4, verbose=FALSE)
  n_rows = base::nrow(iris2d$embed)
  n_cols = base::ncol(iris2d:embed)
  
  expect_equal(n_rows, base::nrow(dat))
  expect_equal(n_cols, 2)
})
