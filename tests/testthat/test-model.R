test_that("model run agrees with reference", {
  pars <- reference_pars()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- c(1, 8, 22, 23)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  reference <- array(
    c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857, 
      0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857, 
      0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857, 
      0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.403355076072801, 0.403355076072801, 
      1, 2.1478536987478, NaN, 0, 0, 8, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.521809241106243, 
      1, 1.93354504572728, NaN, 0, 0, 8, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.489522393762213, 
      1, 2.14825393465371, NaN, 0, 0, 8, 1, 2.21535879122106e-05, 3.3281764902013e-05, 
      0.999945242777042, 3.85493558809803e-06, 1.17221335663037e-05, 
      0.999990350213984, 3.98279297059406e-05, 3.84839746139132e-05, 
      0.999938942777973, 7.98409202313655e-06, 1.2906547416643e-05, 
      0.997545509076154, 4.98560186046243e-06, 1.16454450795262e-06, 
      749958.003941533, 24.9613236765098, 8.7916001747278, 8.24313461576502, 
      2.99456895909997, 21.2216081003379, 2.60318101964684, 0, 5, 0, 
      36.9960584670026, 3.53912804730906, -0.387010010305226, 1, 2.20449400165164, 
      9.18315345100951, 0, 0, 22, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 
      0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.87974458719621, 
      1, 2.16280822646685, NaN, 0, 0, 22, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.13010279522456, 
      1, 2.77860738326335, NaN, 0, 0, 22, 0.999997014307458, 1.95109819258826e-05, 
      3.17954319861917e-05, 0.999942842821874, 3.93729744857565e-06, 
      1.29518207874794e-05, 0.999990350213984, 3.92199679055253e-05, 
      3.76727850380659e-05, 0.999935144523607, 8.63040118619906e-06, 
      1.44904073892214e-05, 0.997413014570077, 5.69548198395564e-06, 
      1.49660861352227e-06, 749955.004419648, 24.8406801016762, 9.71386559060957, 
      10.441034659447, 6.11473441866369, 20.5684289725258, 5.25588203218908, 
      1, 5, 0, 38.9955803517327, 3.15211803700383, -0.387010010305226, 
      1, 2.20449400165164, 8.22478166928513, 0, 0, 23, 1, 0, 0, 1, 
      0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, -1.87974458719621, 1, 2.16280822646685, NaN, 0, 0, 
      23, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, -2.13010279522456, 1, 2.77860738326335, 
      NaN, 0, 0, 23), dim = c(34L, 3L, 4L),
    dimnames = list(
      c("thetaf", "MSEf", "MEf", "MSSf", "MSIf", "MIf", "thetag",
        "MSEg", "MEg", "MSSg", "MSIg", "MIg", "thetah", "MEh", "MIh",
        "S", "E", "I", "R", "newI", "Eseed", "newIseed", "cutf",
        "cutg", "cuth", "cuts", "seedrate", "dseedrate", "theta_vacc",
        "beta", "cumulative_partners", "V1", "V2", "time"), NULL, NULL)
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
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0)

})

test_that("no infections when seedrate = 0", {
  pars <- reference_pars()
  pars$seedrate0 <- pars$seedrate_sd <- 0

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  expect_true(all(res["newI", , ] == 0))
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0)
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
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0)
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
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0)
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
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0)
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


test_that("vaccination works", {
  pars <- reference_pars()
  pars$beta0 <- 2
  pars$seedrate0 <- 2
  n_par <- 3
  t <- seq(1, 365)

  # No vaccination when doses = 0
  pars$vacc_doses <- pars$vacc_doses2 <- 0
  pars$vacc_start_day <- pars$vacc_start_day2 <- 0
  m <- model$new(pars, 1, n_par, seed = 1)
  res_novax <- m$simulate(t)
  rownames(res_novax) <- names(m$info()$index)
  expect_true(all(res_novax["V1", , ] == 0))
  expect_true(all(res_novax["V2", , ] == 0))
  
  # 100% vaccination on day 1 with 99% efficacy
  pars$vacc_doses <- pars$vacc_doses2 <- pars$N
  pars$vacc_start_day <- 1
  pars$vacc_duration <- 1
  pars$vacc_start_day2 <- 2
  pars$vacc_duration2 <- 2
  pars$vacc_efficacy <- pars$vacc_efficacy2 <- 0.99

  m <- model$new(pars, time = 1, n_par, seed = 1)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)
  expect_true(all(res["V1", , 1] == 0))
  expect_true(all(res["V1", , -1] == 1))
  expect_true(all(res["V2", , 1:2] == 0))
  expect_true(all(res["V2", , 3] == 0.5))
  expect_true(all(res["V2", , -(1:3)] == 1))
  
  ## random vaccination
  # pars$vacc_targetted <- 0
  m <- model$new(pars, time = 1, n_par, seed = 1)
  res_random <- m$simulate(t)
  rownames(res_random) <- names(m$info()$index)
  expect_true(all(res_random["V1", ,-1] == 1))
  expect_true(all(res_random["V2",1 ,-(1:3)] == 1))
  
  ## perfectly targeted vaccination
  pars$vacc_targetted <- 1
  m <- model$new(pars, time = 1, n_par, seed = 1)
  res_perfect <- m$simulate(t)
  rownames(res_perfect) <- names(m$info()$index)
  expect_true(all(res_perfect["V1", ,-1] == 1))
  expect_true(all(res_perfect["V2",1 ,-(1:3)] == 1))

  # targetting has no effect when everyone is vaccinated
  expect_equal(res_perfect, res_random)
  
  ## 100% doses day 1, efficacy 0%
  pars$vacc_efficacy <- pars$vacc_efficacy2 <- 0
  m <- model$new(pars, time = 1, n_par, seed = 1)
  res_ve0 <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_ve0) <- names(m$info()$index)
  expect_true(all(res_ve0["V1", ,-1] == 1))
  expect_true(all(res_ve0["V2",1 ,-(1:3)] == 1))
  
  # Check we get the same result when efficacy is 0 vs when there is no vaccination
  expect_equal(res_novax["newI", ,], res_ve0["newI", , ])
  
  pars$vacc_doses <- pars$vacc_doses2 <- 1e5
  pars$vacc_targetted <- 0.8
  pars$vacc_efficacy <- pars$vacc_efficacy2 <- 0.99
  m <- model$new(pars, time = 1, n_par, seed = 1)
  res_V50k <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_V50k) <- names(m$info()$index)
  
  pars$vacc_targetted <- 1
  m <- model$new(pars, time = 1, n_par, seed = 1)
  res_V50k_perfect <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_V50k_perfect) <- names(m$info()$index)
  
  pars$vacc_targetted <- 0
  m <- model$new(pars, time = 1, n_par, seed = 1)
  res_V50k_random <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_V50k_random) <- names(m$info()$index)
  
  res_V50k_perfect - res_V50k_random ## these should not be equal - looks like targetting is not really affecting anything
  #  Plots to investigate
  # legend <- list(ve0 = 2, ve99vt80 = 3, V50kvt0 = 4, V50kvt100 = 1, V50kvt80 = 7)
  # par(mfrow=c(1, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
  # matplot(t(res_novax["newI", ,]), type = "l", col = legend$ve0, lty = 1)
  # matlines(t(res["newI", ,]), col = legend$ve99vt80, lty = 1)
  # matlines(t(res_V50k["newI", ,]), col = legend$V50kvt80, lty = 1)
  # matlines(t(res_V50k_perfect["newI", ,]), col = legend$V50kvt100, lty = 2)
  # matlines(t(res_V50k_random["newI", ,]), col = legend$V50kvt0, lty = 3)
  # legend("topleft", legend = names(legend), fill = unlist(legend))
  
  ## New odin debugging facility
  ## See https://mrc-ide.github.io/odin/articles/debugging.html
  # model_debug <- odin.dust::odin_dust("inst/odin/model.R", debug = TRUE)
  # m <- model_debug$new(pars, time = 1, n_par, seed = 1)
  # res <- m$simulate(t)

})
