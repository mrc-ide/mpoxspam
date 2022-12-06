test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("ll_nbinom behaves in corner cases", {
  expect_equal(ll_nbinom(0, 0, 0.5, Inf), 0)
  expect_equal(ll_nbinom(0, 0, 0, Inf), 0)
  expect_equal(ll_nbinom(0, 0, 1, Inf), 0)
  expect_true(is.finite(ll_nbinom(0, 0, 0.5, 1e6)))
  expect_true(is.finite(ll_nbinom(0, 0, 0, 1e6)))
  expect_true(is.finite(ll_nbinom(0, 0, 1, 1e6)))
})


test_that("ll_binom behaves in corner cases", {
  expect_equal(ll_binom(0, 0, 0, 0, Inf), NaN) # 0 / 0 = undefined
  expect_equal(ll_binom(0, 0, 0, 0, 1e100), 0) # eps / 2eps ~= 0
  expect_equal(ll_binom(0, 0, 1, 1, Inf), 0)
  expect_equal(ll_binom(0, 0, 1, 0, Inf), 0)
  expect_true(is.finite(ll_binom(0, 0, 1, 1, 1e6)))
  expect_true(is.finite(ll_binom(0, 0, 0, 0, 1e6)))
  expect_true(is.finite(ll_binom(0, 0, 1, 0, 1e6)))

  set.seed(1)
  ll1 <- replicate(1000, ll_binom(1, 2, 1, 2, 1e4))
  ll2 <- replicate(1000, ll_binom(1, 2, 1, 2, 1e8))
  expect_lt(var(ll2), var(ll1))
  expect_equal(mean(ll2), mean(ll1), tolerance = 1e-4)
})
