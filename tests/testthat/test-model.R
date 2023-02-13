test_that("model run agrees with reference", {
  pars <- reference_pars()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- c(1, 8, 22, 23)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  reference <- array(
    c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1, 1.97848857973857,
      0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1, 1.97848857973857,
      0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1, 1.97848857973857,
      0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0.403355076072801, 0.403355076072801,
      1, 1, 2.1478536987478, NaN, 8, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.521809241106242,
      1, 1, 1.93354504572728, NaN, 8, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.489522393762213,
      1, 1, 2.14825393465371, NaN, 8, 1, 2.21535879122106e-05, 3.3281764902013e-05,
      0.999945242777042, 3.85493558809803e-06, 1.17221335663037e-05,
      0.999990350213984, 3.98514815773126e-05, 3.85072612911673e-05,
      0.999938916431804, 7.9874185101396e-06, 1.29098740848222e-05,
      0.997545509076154, 4.98560186046243e-06, 1.16454450795262e-06,
      749958, 24.9641510327492, 8.79231888848881, 8.24353007876198,
      2.99487871896417, 21.2244354565773, 2.60349077951105, 0, 5, 0,
      37, 3.53912804730906, -0.387010010305226, 1, 1, 2.20449400165164,
      9.18340590136729, 22, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.87974458719621,
      1, 1, 2.16280822646685, NaN, 22, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.13010279522456,
      1, 1, 2.77860738326335, NaN, 22, 0.999997014307458, 2.02794486924686e-05,
      3.25639408341882e-05, 0.999942074355108, 3.93729448602354e-06,
      1.29518207874794e-05, 0.999990350213984, 3.92408058928984e-05,
      3.76931608806633e-05, 0.999935118441006, 8.63482568625035e-06,
      1.44958132250125e-05, 0.997413014570077, 5.69548198395564e-06,
      1.49660861352227e-06, 749955, 24.8436321536556, 9.71475804546026,
      10.4416098008842, 6.11539759805783, 20.5713810245052, 5.25654521158322,
      1, 5, 0, 39, 3.15211803700383, -0.387010010305226, 1, 1, 2.20449400165164,
      8.22517858000274, 23, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.87974458719621,
      1, 1, 2.16280822646685, NaN, 23, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.13010279522456,
      1, 1, 2.77860738326335, NaN, 23), dim = c(33L, 3L, 4L),
    dimnames = list(
      c("thetaf", "MSEf", "MEf", "MSSf", "MSIf", "MIf", "thetag",
        "MSEg", "MEg", "MSSg", "MSIg", "MIg", "thetah", "MEh", "MIh",
        "S", "E", "I", "R", "newI", "Eseed", "newIseed", "cutf",
        "cutg", "cuth", "cuts", "seedrate", "dseedrate", "theta_vacc",
        "S_vacc", "beta", "cumulative_partners", "time"), NULL, NULL)
    )

  expect_equal(res, reference)

})

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
