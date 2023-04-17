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

  expect_equal(test_h(1, hshape, hrate), evaluate_pgf(h, 1))
  expect_equal(test_h(0.342, hshape, hrate), evaluate_pgf(h, 0.342))
})

test_that("can compute theta_vacc", {
  hshape <- 0.26
  hrate <- 1.85 * 7
  res1 <- reference_update_theta_vacc_4.2(1, 0.00063)
  expect_equal(res1[[1]], res1[[2]], tolerance = 1e-5)
  expect_equal(test_update_theta_vacc4_2(1, 0.00063, hshape, hrate), res1[[2]], tolerance = 1e-5)

  res2 <- reference_update_theta_vacc_4.2(0.96950, 0.00063)
  expect_equal(res2[[1]], res2[[2]], tolerance = 1e-4)
  expect_equal(test_update_theta_vacc4_2(0.96950, 0.00063, hshape, hrate), res2[[2]], tolerance = 1e-5)

  res1 <- reference_update_theta_vacc_4.3(0.5)
  expect_equal(res1[[1]], res1[[2]], tolerance = 1e-5)
  expect_equal(test_update_theta_vacc4_3(0.5, hshape, hrate), res1[[2]], tolerance = 1e-5)
})


test_that("ll_nbinom behaves in corner cases", {
  expect_equal(ll_nbinom(0, 0, 0.5, Inf), 0)
  expect_equal(ll_nbinom(0, 0, 0, Inf), 0)
  expect_equal(ll_nbinom(0, 0, 1, Inf), 0)
  expect_true(is.finite(ll_nbinom(0, 0, 0.5, 1e6)))
  expect_true(is.finite(ll_nbinom(0, 0, 0, 1e6)))
  expect_true(is.finite(ll_nbinom(0, 0, 1, 1e6)))

  expect_equal(ll_nbinom(NA, 0, 0.5, Inf), 0)
})


test_that("ll_betabinom behaves in corner cases", {
  rho <- 0.5
  big <- 1e100
  expect_equal(ll_betabinom(0, 0, 0, 0, rho, Inf), 0)
  expect_equal(ll_betabinom(0, 0, 0, 0, rho, big), 0)
  expect_equal(ll_betabinom(0, 0, 1, 1, rho, Inf), 0)

  ## These could probably be special cases, they fail because of Inf -
  ## Inf in the beta fraction
  expect_equal(ll_betabinom(0, 0, 1, 0, rho, Inf), 0)
  expect_equal(ll_betabinom(0, 0, 0, 1, rho, Inf), 0)
  expect_equal(ll_betabinom(0, 0, 1, 1, rho, Inf), 0)
  expect_equal(ll_betabinom(0, 0, 1, 0, rho, big), 0)
  expect_equal(ll_betabinom(0, 0, 0, 1, rho, big), 0)
  expect_equal(ll_betabinom(0, 0, 1, 1, rho, big), 0)

  set.seed(1)
  ll1 <- replicate(1000, ll_betabinom(1, 2, 1, 2, rho, 1e4))
  ll2 <- replicate(1000, ll_betabinom(1, 2, 1, 2, rho, 1e8))
  expect_lt(var(ll2), var(ll1))
  expect_equal(mean(ll2), mean(ll1), tolerance = 1e-4)
})


test_that("ll_nbinom behaves in corner cases", {
  rng <- dust:::dust_rng_pointer$new(seed = 1L)
  expect_equal(test_ll_nbinom(0, 0, 0.5, Inf, rng), 0)
  expect_equal(test_ll_nbinom(0, 0, 0, Inf, rng), 0)
  expect_equal(test_ll_nbinom(0, 0, 1, Inf, rng), 0)
  expect_true(is.finite(test_ll_nbinom(0, 0, 0.5, 1e6, rng)))
  expect_true(is.finite(test_ll_nbinom(0, 0, 0, 1e6, rng)))
  expect_true(is.finite(test_ll_nbinom(0, 0, 1, 1e6, rng)))

  expect_equal(test_ll_nbinom(NA, 0, 0.5, Inf, rng), 0)
})


test_that("compiled ll_betabinom agrees", {
  rng <- dust:::dust_rng_pointer$new(seed = 1L)

  rho <- 0.5
  big <- 1e100
  expect_equal(test_ll_betabinom(0, 0, 0, 0, rho, Inf, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 0, 0, rho, big, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 1, 1, rho, Inf, rng), 0)

  ## These could probably be special cases, they fail because of Inf -
  ## Inf in the beta fraction
  expect_equal(test_ll_betabinom(0, 0, 1, 0, rho, Inf, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 0, 1, rho, Inf, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 1, 1, rho, Inf, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 1, 0, rho, big, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 0, 1, rho, big, rng), 0)
  expect_equal(test_ll_betabinom(0, 0, 1, 1, rho, big, rng), 0)

  ll1 <- replicate(1000, test_ll_betabinom(1, 2, 1, 2, rho, 1e4, rng))
  ll2 <- replicate(1000, test_ll_betabinom(1, 2, 1, 2, rho, 1e8, rng))
  expect_lt(var(ll2), var(ll1))
  expect_equal(mean(ll2), mean(ll1), tolerance = 1e-4)

  expect_equal(ll_betabinom(2, 4, 89, 129, 0.5, Inf),
               test_ll_betabinom(2, 4, 89, 129, 0.5, Inf, rng))
})


test_that("hu* functions behave with very small thetav", {
 
  test_hu_theta <- function(thetav) {
    x <- 2.7582643254561925e-05
    vr <- 0.80000000000000004
    V1 <- 0.30823038793678742
    V2 <- 0.0080929608005416483
    v1eff <- 0.89000000000000001
    v2eff <- 0.94999999999999996
    hs <- 0.066370290188591391
    hr <- 3.3057509920856094
    test_hu(x, vr, V1, V2, v1eff, v2eff, thetav, hs, hr)
  }
  thetav <- 4.0732983338996529e-157
  thetav <- exp(seq(log(1e-2), log(1e-150), length.out = 100))
  y <- sapply(thetav, test_hu_theta)
  
  par(mfrow = c(2, 2))
  for (i in 1:4) {
  plot(thetav, y[i, ], log = "x", type = "l",
       ylab = sprintf("hu%s(thetav)", c("", "p", "pp", "ppp")[i]))
  }

  eps <- .Machine$double.eps
  expect_equal(test_hu_theta(eps ^ (1/3)) == test_hu_theta(eps), c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(test_hu_theta(eps ^ (1/2)) == test_hu_theta(eps), c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(test_hu_theta(eps* 1.1) == test_hu_theta(eps), c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(test_hu_theta(eps) == test_hu_theta(eps), c(TRUE, TRUE, TRUE, TRUE))
  expect_false(any(test_hu_theta(eps ^ (1/3) * 2) == test_hu_theta(eps ^ (1/3))))

})
