test_that("tSNE", {
  myn = 10
  myp = 3
  
  prep_mat = rbind(matrix(rnorm(myn*myp, mean=-2), ncol=myp),
                   matrix(rnorm(myn*myp, mean= 2), ncol=myp))
  run_tsne = maotai::tsne(prep_mat, ndim=2, perplexity=5)$embed

  expect_equal(base::nrow(run_tsne), 2*myn)
  expect_equal(base::ncol(run_tsne), 2)
})
