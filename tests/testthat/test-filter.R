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
  expect_equal(ll, -453.653978547233)
})
