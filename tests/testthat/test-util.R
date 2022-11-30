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
  expect_equal(ll_binom(0, 0, 0, Inf), 0)
  expect_equal(ll_binom(0, 0, 0.5, Inf), 0)
  expect_equal(ll_binom(0, 0, 1, Inf), 0)
  expect_true(is.finite(ll_binom(0, 0, 0.5, 1e6)))
  expect_true(is.finite(ll_binom(0, 0, 0, 1e6)))
  # expect_true(is.finite(ll_binom(0, 0, 1, 1e6))) # error

  expect_true(is.finite(ll_binom(1, 2, 0, 1e6)))
})
