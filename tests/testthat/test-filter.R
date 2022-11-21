test_that("can run filter", {
  pars <- reference_pars()
  dat <- filter_data()
  filter <- m4_2_filter(dat, n_particles = 5000, n_threads = 1, seed = 1L)
  ## pomp:   :  207s
  ## 8 threads: 0.896 user, 0.114 elapsed (1725x pomp)
  ## 1 thread:  0.665 user, 0.664 elapsed ( 259x pomp)
  system.time(filter$run(pars))
})