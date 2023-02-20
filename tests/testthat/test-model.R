test_that("model run agrees with reference", {
  pars <- reference_pars()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- c(1, 8, 22, 23)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  reference <- array(
    c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857,
      0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857,
      0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857,
      0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.403355076072801, 0.403355076072801,
      1, 2.1478536987478, NaN, 0, 0, 8, 1, 0, 0, 1, 0, 0, 1, 0, 0,
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.521809241106242,
      1, 1.93354504572728, NaN, 0, 0, 8, 1, 0, 0, 1, 0, 0, 1, 0, 0,
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.489522393762213,
      1, 2.14825393465371, NaN, 0, 0, 8, 1, 2.21564447569044e-05, 3.32855347103323e-05,
      0.999945237595088, 3.85536371982324e-06, 1.17230918513184e-05,
      0.999990350213984, 3.98553438687456e-05, 3.85110310750159e-05,
      0.999938911321708, 7.98807244863049e-06, 1.29108323663166e-05,
      0.99754524783798, 4.98560180691493e-06, 1.16454450026421e-06,
      749958, 24.9641510327492, 8.79231888848881, 8.24353007876198,
      2.99487871896417, 21.2244354565773, 2.60349077951105, 0, 5, 0,
      37, 3.53912804730906, -0.387010010305226, 1, 2.20449400165164,
      9.18340578319603, 0, 0, 22, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
      0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.87974458719621,
      1, 2.16280822646685, NaN, 0, 0, 22, 1, 0, 0, 1, 0, 0, 1, 0, 0,
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.13010279522456,
      1, 2.77860738326335, NaN, 0, 0, 22, 0.999997014307458, 1.95137854115229e-05,
      3.17993680554976e-05, 0.999942837028392, 3.93778084572668e-06,
      1.29530107272803e-05, 0.999990350213984, 3.92448596010161e-05,
      3.76970969197487e-05, 0.999935112744579, 8.63560560133797e-06,
      1.44970031591145e-05, 0.997412721681836, 5.69548191014851e-06,
      1.49660860106252e-06, 749955, 24.8436321536556, 9.71475804546026,
      10.4416098008842, 6.11539759805783, 20.5713810245052, 5.25654521158322,
      1, 5, 0, 39, 3.15211803700383, -0.387010010305226, 1, 2.20449400165164,
      8.22517844981633, 0, 0, 23, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
      0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.87974458719621,
      1, 2.16280822646685, NaN, 0, 0, 23, 1, 0, 0, 1, 0, 0, 1, 0, 0,
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.13010279522456,
      1, 2.77860738326335, NaN, 0, 0, 23), dim = c(34L, 3L, 4L),
    dimnames = list(
      c("thetaf", "MSEf", "MEf", "MSSf", "MSIf", "MIf", "thetag",
        "MSEg", "MEg", "MSSg", "MSIg", "MIg", "thetah", "MEh", "MIh",
        "S", "E", "I", "R", "newI", "Eseed", "newIseed", "cutf",
        "cutg", "cuth", "cuts", "seedrate", "dseedrate", "theta_vacc",
        "beta", "cumulative_partners", "V1", "V2", "time"), NULL, NULL)
    )

  expect_equal(res, reference)

}

test_that("when beta = 0 all infections are from travel", {
  pars <- reference_pars()
  pars$beta0 <- pars$beta_sd <- 0
  pars$seedrate0 <- 10

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_equal(res["newI", , ], res["newIseed", , ])
  expect_true(all(colSums(res[c("S", "E", "I", "R"), , ]) == pars$N))

})

test_that("no infections when seedrate = 0", {
  pars <- reference_pars()
  pars$seedrate0 <- pars$seedrate_sd <- 0

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_true(all(res["newI", , ] == 0))
  expect_true(all(colSums(res[c("S", "E", "I", "R"), , ]) == pars$N))
})

test_that("infections but no seeding when seedrate = 0, i0 > 0", {
  pars <- reference_pars()
  pars$seedrate0 <- pars$seedrate_sd <- 0
  pars$i0 <- 10

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_true(all(res["newIseed", , ] == 0))
  expect_true(any(res["newI", , ] > 0))
  expect_true(all(colSums(res[c("S", "E", "I", "R"), , ]) == pars$N))
})

test_that("no-one leaves E when gamma0 = 0", {
  pars <- reference_pars()
  pars$gamma0 <- 0
  pars$beta0 <- 100
  pars$seedrate0 <- 100


  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_true(all(res["I", , ] == 0))
  expect_true(all(res["R", , ] == 0))
  expect_true(all(res["E", , -1] > 0))
  expect_true(all(apply(res["E", , ], 1, diff) > 0))
  expect_true(all(colSums(res[c("S", "E", "I", "R"), , ]) == pars$N))
})

test_that("no-one leaves I when gamma1 = 0", {
  pars <- reference_pars()
  pars$gamma1 <- 0
  pars$i0 <- 1e5
  pars$beta0 <- pars$beta_sd <- 0
  pars$seedrate0 <- pars$seedrate_sd <- 0


  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_true(all(res["R", , ] == 0))
  expect_true(all(res["I", , -1] > 0))
  expect_true(all(apply(res["I", , ], 1, diff) > 0))
  expect_true(all(apply(res["E", , ], 1, diff) < 0))
  expect_true(all(colSums(res[c("S", "E", "I", "R"), , ]) == pars$N))
})

test_that("user input fixed beta works", {
  pars <- reference_pars()
  t <- seq(1, 21)

  pars$beta_step[] <- pars$beta0
  pars$dseedrate_step[] <- 0
  pars$stochastic_behaviour <- 0

  m_det <- model$new(pars, 1, 3, seed = 1)
  res_det <- m_det$simulate(t)
  rownames(res_det) <- names(m_det$info()$index)

  pars$stochastic_behaviour <- 1
  pars$beta_sd <- pars$seedrate_sd <- 0

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_equal(res, res_det)

  ## Problem: need to specify a [nsample, nstep] array for each time varying
  ## parameter for this to be useful - probably better done in the simulation
  ## task, as we'll need to break it apart every 7 days anyway
  pars$beta_sd <- pars$seedrate_sd <- 2

  m <- model$new(pars, 1, 3, seed = 2)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  # check can input betas and dseedrate for a single particle
  pars$beta_step <- res["beta", 1, ]
  pars$dseedrate_step <- res["dseedrate", 1, ]
  pars$stochastic_behaviour <- 0

  m_det <- model$new(pars, 1, 3, seed = 2)
  t <- seq(1, 21)
  res_det <- m_det$simulate(t)
  rownames(res_det) <- names(m_det$info()$index)

  expect_equal(res_det[, 1, ], res[, 1, ])
  expect_equal(res_det[, -1, 1], res[, -1, 1])
  expect_false(all(res_det[, -1, -1] == res[, -1, -1]))

})
