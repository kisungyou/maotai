library(maotai)

test_that("cov2corr works", {
  myp = 5
  myn = 50
  
  prep_mat = stats::cov(matrix(rnorm(myn*myp), ncol=myp))
  vec_corr = diag(maotai::cov2corr(prep_mat))
  
  expect_equal(vec_corr, rep(1, myp))
})

test_that("cov2pcorr works", {
  myp = 5
  myn = 50
  
  prep_mat = stats::cov(matrix(rnorm(myn*myp), ncol=myp))
  vec_pcor = diag(maotai::cov2pcorr(prep_mat))
  
  expect_equal(vec_pcor, rep(1, myp))
})