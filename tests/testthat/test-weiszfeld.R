test_that("weiszfeld", {
  n = 20
  p = 10
  t = seq(from=0,to=10,length.out=p)
  
  X = array(0,c(n,p))
  for (i in 1:n){
     X[i,] = sin(t) + stats::rnorm(p, sd=0.5)
  }

  vecL1 = as.vector(weiszfeld(X))
  expect_equal(length(vecL1), p)
})
