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

