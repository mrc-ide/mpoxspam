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
  pars$compare_cases <- "negbinom"
  pars$compare_travel <- "betabinom"
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

test_that("error on unreognised observation distribution", {

  dat <- reference_data()
  pars <- reference_pars()
  pars$compare_travel <- "foo"

  m <- model$new(pars, 1, 100, seed = 1)
  m$set_index(model_index(m$info())$run)
  y <- m$run(dat$day[[2]])

  expect_error(model_compare(y, dat[2, ], pars),
               "unrecognised compare function foo")

  pars$compare_cases <- "bar"

  m <- model$new(pars, 1, 100, seed = 1)
  m$set_index(model_index(m$info())$run)
  y <- m$run(dat$day[[2]])

  expect_error(model_compare(y, dat[2, ], pars),
               "unrecognised compare function bar")
})

test_that("can use negbinom dist for travel", {
  dat <- reference_data()
  pars <- reference_pars()
  pars$kappa_travel <- pars$kappa_cases
  pars$compare_travel <- "negbinom"

  m <- model$new(pars, 1, 100, seed = 1)
  m$set_index(model_index(m$info())$run)

  y <- m$run(dat$day[[2]])
  expect_equal(sum(model_compare(y, dat[2, ], pars)), -984.83445)
})
