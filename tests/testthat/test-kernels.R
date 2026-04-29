test_that("uniform kernel gives 1 inside bandwidth, 0 outside", {
  x   <- c(-0.6, -0.4, 0, 0.4, 0.6)
  bw  <- 0.5
  w   <- kernel_weights(x, bw, "uniform")
  expect_equal(w, c(0, 1, 1, 1, 0))
})

test_that("triangular kernel at boundary is 0, at center is bwidth", {
  bw  <- 0.5
  w   <- kernel_weights(c(-bw, 0, bw), bw, "triangular")
  expect_equal(w[1], 0)
  expect_equal(w[2], bw)
  expect_equal(w[3], 0)
})

test_that("epanechnikov kernel is non-negative and zero at boundary", {
  bw  <- 0.5
  x   <- c(-0.5, -0.25, 0, 0.25, 0.5)
  w   <- kernel_weights(x, bw, "epanechnikov")
  expect_true(all(w >= 0))
  expect_equal(w[c(1, 5)], c(0, 0))
})

test_that("obs_weights multiply into kernel weights", {
  x  <- c(-0.3, 0, 0.3)
  ow <- c(2, 3, 4)
  w_base <- kernel_weights(x, 0.5, "uniform")
  w_mult <- kernel_weights(x, 0.5, "uniform", obs_weights = ow)
  expect_equal(w_mult, w_base * ow)
})
