library(maotai)

test_that("Bayesian MDS", {
  myn = 10
  myp = 3
  
  prep_mat = rbind(matrix(rnorm(myn*myp, mean=-2), ncol=myp),
                   matrix(rnorm(myn*myp, mean= 2), ncol=myp))
  
  iris2d = maotai::bmds(prep_mat, ndim=2, mc.iter=2, par.step=4, verbose=FALSE)
  n_rows = base::nrow(iris2d$embed)
  n_cols = base::ncol(iris2d$embed)
  
  expect_equal(n_rows, base::nrow(prep_mat))
  expect_equal(n_cols, 2)
})
