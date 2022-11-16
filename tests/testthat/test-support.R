test_that("f functions agree", {
  f <- list(function(x) (2809 + 2000 * x + 95 * x^2) / 4904,
            function(x) (2000 + 2 * 95 * x) / 4904,
            function(x) (2 * 95) / 4904,
            function(x) 0)

  expect_equal(test_f(1), evaluate_pgf(f, 1))
  expect_equal(test_f(0.342), evaluate_pgf(f, 0.342))
})

test_that("g functions agree", {
  g <- list(function(x) (2943 + 1009 * x + 477 * x^2 + 475 * x^3) / 4904,
            function(x) (1009 + 2 * 477 * x^1 + 3 * 475 * x^2) / 4904,
            function(x) (2 * 477 + 2 * 3 * 475 * x^1) / 4904,
            function(x) (2 * 3 * 475) / 4904)

  expect_equal(test_g(1), evaluate_pgf(g, 1))
  expect_equal(test_g(0.342), evaluate_pgf(g, 0.342))
})

test_that("h functions agree", {
  hshape <- 0.26
  hrate <- 1.85 * 7

  h <- list(function(x) (1 - log(x) / hrate)^(-hshape),
            function(x) hshape * (1 - log(x) / hrate)^(-hshape - 1) / ((hrate * x)),
            function(x) hshape * ((hrate - log(x)) / hrate)^(-hshape) * (-hrate + hshape + log(x) + 1) / (x^2 * (hrate - log(x))^2),
            function(x) hshape * ((hrate - log(x)) / hrate)^(-hshape) * (hshape^2 + 3 * hshape + 2 * (hrate - log(x))^2 - 3 * (hrate - log(x)) * (hshape + 1) + 2) / (x^3 * (hrate - log(x))^3))

  expect_equal(test_h(1), evaluate_pgf(h, 1))
  expect_equal(test_h(0.342), evaluate_pgf(h, 0.342))
})
