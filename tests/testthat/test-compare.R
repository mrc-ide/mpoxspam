test_that("compare agrees with reference", {
  pars <- reference_pars()
  dat <- reference_data()

  m <- m4_2$new(pars, 1, 100, seed = 1)
  m$set_index(m4_2_index(m$info())$run)

  for (i in seq_len(nrow(dat))) {
    y <- m$run(dat$day[[i]])
    expect_equal(m4_2_compare(y, dat[i, ], pars), reference_compare(y, dat[i, ], pars))
  }
})


test_that("can use compiled compare", {
  pars <- reference_pars()
  dat <- filter_data()

  m <- m4_2$new(pars, 1, 100, seed = 1)
  m$set_index(m4_2_index(m$info())$run)
  m$set_data(dust::dust_data(dat, "time_end"))

  for (i in seq_len(nrow(dat))) {
    y <- m$run(dat$day_end[[i]])
    expect_equal(m$compare_data(),
                 m4_2_compare(y, dat[i, ], pars))
  }
})
