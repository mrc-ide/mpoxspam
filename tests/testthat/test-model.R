test_that("model run agrees with reference", {
  pars <- reference_pars()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- c(1, 8, 22, 23)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  reference <- array(
    c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857, 
      0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 
      0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.111546627766835, 
      0, 1, 1.97848857973857, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 
      1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, -0.111546627766835, 0, 1, 1.97848857973857, 0, 0, 
      0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 
      0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.97848857973857, 
      NaN, 0, 0, NaN, NaN, NaN, NaN, 1, 1, 0.8, 1, 0, 0, 1, 0, 0, 1, 
      0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 1, 1.97848857973857, NaN, 0, 0, NaN, NaN, NaN, NaN, 1, 
      1, 0.8, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.97848857973857, NaN, 
      0, 0, NaN, NaN, NaN, NaN, 1, 1, 0.8, 1, 0, 0, 1, 0, 0, 1, 0, 
      0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 1, 1.97848857973857, NaN, 0, 0, NaN, NaN, NaN, NaN, 1, 1, 
      2.2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.97848857973857, NaN, 0, 
      0, NaN, NaN, NaN, NaN, 1, 1, 2.2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
      1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      1, 1.97848857973857, NaN, 0, 0, NaN, NaN, NaN, NaN, 1, 1, 2.2, 
      1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.97848857973857, NaN, 0, 0, NaN, 
      NaN, NaN, NaN, 1, 1, 2.3, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
      1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.97848857973857, 
      NaN, 0, 0, NaN, NaN, NaN, NaN, 1, 1, 2.3, 1, 0, 0, 1, 0, 0, 1, 
      0, 0, 1, 0, 0, 1, 0, 0, 750000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 1, 1.97848857973857, NaN, 0, 0, NaN, NaN, NaN, NaN, 1, 
      1, 2.3), dim = c(40L, 3L, 4L),
    dimnames = list(c("thetaf", "MSEf", 
     "MEf", "MSSf", "MSIf", "MIf", "thetag", "MSEg", "MEg", "MSSg", 
     "MSIg", "MIg", "thetah", "MEh", "MIh", "S", "E", "I", "R", "newI", 
     "Eseed", "newIseed", "cutf", "cutg", "cuth", "cuts", "seedrate", 
     "dseedrate", "theta_vacc", "beta", "cumulative_partners", "V1", 
     "V2", "Reff_f", "Reff_g", "Reff_h", "Reff", "vacc_period", "vacc_period2", 
     "time"), NULL, NULL))

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
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0)

})


test_that("vaccination works", {
  pars <- reference_pars()
  pars$beta0 <- 2
  pars$seedrate0 <- 2
  n_par <- 3
  t <- seq(10, 14 * 10, by = 10)

  # No vaccination when doses = 0
  pars$vacc_doses_step <- pars$vacc_doses2_step <- 0
  pars$vacc_start_day_step <- pars$vacc_start_day2_step <- 0
  m <- model$new(pars, 0, n_par, seed = 1)
  res_novax <- m$simulate(t)
  rownames(res_novax) <- names(m$info()$index)
  expect_true(all(res_novax["V1", , ] == 0))
  expect_true(all(res_novax["V2", , ] == 0))
  expect_equal(sum(colSums(res_novax[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  
  # 80% vaccination on day 1 with 99% efficacy
  pars$vacc_doses_step <- pars$vacc_doses2_step <- pars$N 
  pars$vacc_start_day_step <- 1
  pars$vacc_duration_step <- 1
  pars$vacc_start_day2_step <- 2
  pars$vacc_duration2_step <- 2
  pars$vacc_efficacy <- pars$vacc_efficacy2 <- 0.80

  m <- model$new(pars, time = 0, n_par, seed = 1)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)
  expect_true(all(res["V1", , 1] == 0))
  expect_equal(sum(res["V1", , -1] - 1), 0, tolerance = 1e-6)
  expect_true(all(res["V2", , 1:2] == 0))
  expect_equal(sum(res["V2", , 3] - 0.5), 0, tolerance = 1e-6)
  expect_equal(sum(res["V2", , -(1:3)] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  
  ## random vaccination
  pars$vacc_targetted <- 0
  m <- model$new(pars, time = 0, n_par, seed = 1)
  res_random <- m$simulate(t)
  rownames(res_random) <- names(m$info()$index)
  expect_equal(sum(res_random["V1", , -1] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(res_random["V2", , -(1:3)] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(colSums(res_random[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  
  ## perfectly targeted vaccination
  pars$vacc_targetted <- 0.8
  m <- model$new(pars, time = 0, n_par, seed = 1)
  res_perfect <- m$simulate(t)
  rownames(res_perfect) <- names(m$info()$index)
  expect_equal(sum(res_perfect["V1", , -1] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(res_perfect["V2", , -(1:3)] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(colSums(res_perfect[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)

  ## 100% doses day 1, efficacy 0%
  pars$vacc_efficacy <- pars$vacc_efficacy2 <- 0
  m <- model$new(pars, time = 0, n_par, seed = 1)
  res_ve0 <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_ve0) <- names(m$info()$index)
  expect_equal(sum(res_ve0["V1", , -1] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(res_ve0["V2", , -(1:3)] - 1), 0, tolerance = 1e-6)
  expect_equal(sum(colSums(res_ve0[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  
  # Check we get the same result when efficacy is 0 vs when there is no vaccination
  expect_equal(res_novax["newI", ,], res_ve0["newI", , ])
  
  pars$vacc_doses_step <- pars$vacc_doses2_step <- 1e5
  pars$vacc_targetted <- 0.8
  pars$vacc_efficacy <- pars$vacc_efficacy2 <- 0.99
  m <- model$new(pars, time = 0, n_par, seed = 1)
  res_V50k <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_V50k) <- names(m$info()$index)
  expect_equal(sum(colSums(res_V50k[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  
  pars$vacc_targetted <- 1
  m <- model$new(pars, time = 0, n_par, seed = 1)
  res_V50k_perfect <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_V50k_perfect) <- names(m$info()$index)
  
  pars$vacc_targetted <- 0
  m <- model$new(pars, time = 0, n_par, seed = 1)
  res_V50k_random <- m$simulate(t)
  ## vaccine efficacy does not affect results
  rownames(res_V50k_random) <- names(m$info()$index)
  
  #  # Plots to investigate
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


test_that("vaccination works after large epidemic",{
  pars <- reference_pars()
  pars$beta0 <- 20
  pars$seedrate0 <- 2
  n_par <- 3
  t <- seq(1, 200)
  
  # No vaccination
  pars$vacc_doses_step <- pars$vacc_doses2_step <- 0
  pars$vacc_start_day_step <- pars$vacc_start_day2_step <- 0
  m <- model$new(pars, 1, n_par, seed = 1)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)
  
  ## Infections still happen when S = 0
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  

  # par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
  # matplot(t(res["R", , ]), type = "l", lty = 1, col = 4, ylab = "S, I, R")
  # matlines(t(res["E", , ]), lty = 1, col = 2)
  # matlines(t(res["I", , ]), lty = 1, col = 3)
  # matlines(t(res["S", , ]), lty = 1, col = 1)
  # matlines(t(res["V1", , ]), lty = 1, col = 5)
  # matlines(t(res["V2", , ]), lty = 1, col = 5)
  # legend("bottomright", legend = c("S", "E", "I", "R", "V1+V2"), fill = 1:5)
  # 
  # matplot(t(colSums(res[c("R", "E", "I", "S"), , ])), type = "l", lty = 1,
  #         col = 1, ylim = c(0, 15e5), ylab = "S+E+I+R")
  # matlines(t(res["E", , ]), lty = 1, col = 2)
  # matlines(t(res["I", , ]), lty = 1, col = 3)
  
  pars <- reference_pars()
  pars$beta0 <- 20
  pars$seedrate0 <- 0
  pars$seedrate_sd <- 0
  pars$beta_sd <- 0
  pars$i0 <- 10
  n_par <- 3
  t <- seq(1, 200)
  
  # No vaccination or seeding
  pars$vacc_doses_step <- pars$vacc_doses2_step <- 0
  pars$vacc_start_day_step <- pars$vacc_start_day2_step <- 0
  m <- model$new(pars, 1, n_par, seed = 1)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)
  
  ## Infections still happen when S = 0
  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)

  
  # par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
  # matplot(t(res["R", , ]), type = "l", lty = 1, col = 4, ylab = "S, I, R")
  # matlines(t(res["E", , ]), lty = 1, col = 2)
  # matlines(t(res["I", , ]), lty = 1, col = 3)
  # matlines(t(res["S", , ]), lty = 1, col = 1)
  # matlines(t(res["V1", , ]), lty = 1, col = 5)
  # matlines(t(res["V2", , ]), lty = 1, col = 5)
  # legend("bottomright", legend = c("S", "E", "I", "R", "V1+V2"), fill = 1:5)
  # 
  # matplot(t(colSums(res[c("R", "E", "I", "S"), , ])), type = "l", lty = 1,
  #         col = 6, ylim = c(0, 15e5), ylab = "S+E+I+R")
  # matlines(t(res["S", , ]), lty = 1, col = 1)
  # 
  # matplot(t(res["S", , ]), type = "l", lty = 1, col = 1, ylab = "S, newI")
  # matlines(t(res["newI", , ]), lty = 1, col = 3)

  
  ## what about with a smaller beta
  
  pars <- reference_pars()
  pars$beta0 <- 5
  pars$beta_sd <- 0
  pars$i0 <- 100
  pars$seedrate0 <- 10
  pars$seedrate_sd <- 1
  t <- seq(1, 800)
  
  # No vaccination
  pars$vacc_doses_step <- pars$vacc_doses2_step <- 0
  pars$vacc_start_day_step <- pars$vacc_start_day2_step <- 0
  m <- model$new(pars, 1, n_par, seed = 1)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)
  
  # par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
  # matplot(t(res["R", , ]), type = "l", lty = 1, col = 4, ylab = "S, I, R")
  # matlines(t(res["E", , ]), lty = 1, col = 2)
  # matlines(t(res["I", , ]), lty = 1, col = 3)
  # matlines(t(res["S", , ]), lty = 1, col = 1)
  # matlines(t(res["V1", , ]), lty = 1, col = 5)
  # matlines(t(res["V2", , ]), lty = 1, col = 5)
  # legend("topright", legend = c("S", "E", "I", "R", "V1+V2"), fill = 1:5)
  # 
  # matplot(t(colSums(res[c("R", "E", "I", "S"), , ])), type = "l", lty = 1,
  #         col = 6, ylim = c(0, 15e5), ylab = "S+E+I+R")
  # matlines(t(res["S", , ]), lty = 1, col = 1)
  # 
  # matplot(t(res["newI", , ]), type = "l", lty = 1, col = 3, ylab = "newI")
  

})

test_that("can implement time-varying vaccination", {
  pars <- reference_pars()
  pars$beta0 <- 2
  pars$seedrate0 <- 2
  n_par <- 2
  t <- seq(0, 4*7 * 10, by = 10)
  
  pars$vacc_doses_step <- c(1e3, 2e3, 0)
  pars$vacc_doses2_step <- c(5e2, 1e3, 0)
  pars$vacc_duration_step <- pars$vacc_duration2_step <- c(7, 7, Inf)
  pars$vacc_start_day_step <- c(0, cumsum(pars$vacc_duration_step) + 1)
  pars$vacc_start_day2_step <- 7 + pars$vacc_start_day_step 


  m <- model$new(pars, 0, n_par, seed = 1)
  res <- m$simulate(t)
  rownames(res) <- names(m$info()$index)

  # par(mfrow = c(1, 2))
  # matplot(t(res["V1", , ]), type = "l", x = t / 10)
  # matplot(t(res["V2", , ]), type = "l", x = t / 10)
  # expect_true(all(res["V1", 1, ] == res["V1", 2, ])) ## doesn't vary by particle

  expect_equal(max(res["V1", 1, ]) * pars$N, sum(pars$vacc_doses_step))
  expect_equal(max(res["V2", 1,]) * pars$N, sum(pars$vacc_doses2_step))
  expect_equal(res["V1", 1, 8] * pars$N, pars$vacc_doses_step[1], ignore_attr = TRUE)
  expect_equal(res["V2", 1, 16] * pars$N, pars$vacc_doses2_step[1], ignore_attr = TRUE)
  expect_equal(res["V2", 1, 23] * pars$N, sum(pars$vacc_doses2_step[1:2]), ignore_attr = TRUE)

  expect_equal(sum(colSums(res[c("S", "E", "I", "R"), , ]) - pars$N), 0, tolerance = 1e-6)
  
})
