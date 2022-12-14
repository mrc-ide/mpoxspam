test_that("can run filter", {
  pars <- reference_pars()
  dat <- filter_data()
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 1, seed = 1L)
  ## pomp:   :  207s
  ## 8 threads: 0.896 user, 0.114 elapsed (1725x pomp)
  ## 1 thread:  0.665 user, 0.664 elapsed ( 259x pomp)
  ll <- filter$run(pars)
  ## Smoke test, will need updating on changes to basically anything,
  ## but guard against unexpected changes.
  expect_equal(ll, -927.94100164507)
})


test_that("can run filter with compiled compare", {
  pars <- reference_pars()
  dat <- filter_data()
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 4, seed = 1L,
                         use_compiled_compare = TRUE)
  ll <- filter$run(pars)
  ## Smoke test, will need updating on changes to basically anything,
  ## but guard against unexpected changes; wildly different to above,
  ## but that's expected with a small number of particles.
  expect_equal(ll, -267.227750645814)
})


test_that("can run filter with negative/beta binomial likelihood", {
  pars <- reference_pars()
  pars$compare_cases <- "negbinom"
  pars$compare_travel <- "betabinom"
  dat <- filter_data()
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 1, seed = 1L)
  ll <- filter$run(pars)
  expect_equal(ll, -274.251863608737)
})

test_that("can run filter with negative/negative binomial likelihood", {
  pars <- reference_pars()
  pars$kappa_travel <- 1 / pars$rho_travel
  pars$compare_cases <- "negbinom"
  pars$compare_travel <- "negbinom"
  dat <- filter_data()
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 1, seed = 1L)
  ll <- filter$run(pars)
  expect_equal(ll, -247.481051609677)
})
