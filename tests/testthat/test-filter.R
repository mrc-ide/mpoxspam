
test_that("can run filter with compiled compare", {
  pars <- reference_pars()
  dat <- filter_data(pars$dt)
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 4, seed = 1L,
                         use_compiled_compare = TRUE)
  ll <- filter$run(pars)
  ## Smoke test, will need updating on changes to basically anything,
  ## but guard against unexpected changes; wildly different to above,
  ## but that's expected with a small number of particles.
  expect_equal(ll, -251.023053997395)
})


test_that("can run filter with negative/beta binomial likelihood", {
  pars <- reference_pars()
  pars$compare_cases <- "negbinom"
  pars$compare_travel <- "betabinom"
  dat <- filter_data(pars$dt)
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 1, seed = 1L)
  ll <- filter$run(pars)
  expect_equal(ll, -269.547690559975)
})

test_that("can run filter with negative/negative binomial likelihood", {
  pars <- reference_pars()
  pars$kappa_travel <- 1 / pars$rho_travel
  pars$compare_cases <- "negbinom"
  pars$compare_travel <- "negbinom"
  dat <- filter_data(pars$dt)
  set.seed(1)
  filter <- model_filter(dat, n_particles = 100, n_threads = 1, seed = 1L)
  ll <- filter$run(pars)
  expect_equal(ll, -251.242574606367)
})

