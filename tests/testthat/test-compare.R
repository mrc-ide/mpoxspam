test_that("compare agrees with reference", {
  pars <- reference_pars()
  dat <- reference_data()

  m <- model$new(pars, 1, 100, seed = 1)
  m$set_index(model_index(m$info())$run)

  for (i in seq_len(nrow(dat))) {
    y <- m$run(dat$day[[i]])
    expect_equal(
      model_compare(y, dat[i, ], pars),
      reference_compare(y, dat[i, ], pars))
  }
})


test_that("compare is correct for missing data", {
  pars <- reference_pars()
  dat <- reference_data()
  dat[2, "Ytravel"] <- NA

  m <- model$new(pars, 1, 100, seed = 1)
  m$set_index(model_index(m$info())$run)
  m$set_data(dust::dust_data(dat, "day"))

  y <- m$run(dat$day[[2]])
  expect_equal(reference_compare(y, dat[2, ], pars), rep(0, 100))
  expect_equal(model_compare(y, dat[2, ], pars), rep(0, 100))
  expect_equal(m$compare_data(), rep(0, 100))
})


test_that("can use compiled compare", {
  pars <- reference_pars()
  pars$use_new_compare <- 1
  pars$exp_noise <- Inf
  dat <- filter_data()

  m <- model$new(pars, 1, 100, seed = 1)
  m$set_index(model_index(m$info())$run)
  m$set_data(dust::dust_data(dat, "time_end"))

  for (i in seq_len(nrow(dat))) {
    y <- m$run(dat$day_end[[i]])
    expect_equal(m$compare_data(),
                 model_compare(y, dat[i, ], pars))
  }
})
