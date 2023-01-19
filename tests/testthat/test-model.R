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
