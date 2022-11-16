compare_m4 <- function(state, observed, pars) {
  newI <- state[, "newI"]
  MSI <- state[, "MSI"]
  MSIg <- state[, "MSIg"]
  MIh <- state[, "MIh"]
  seedrate <- state[, "seedrate"]

  hshape <- 0.26
  hrate <- 1.85 * 7
  fp1 <- (2000 + 2 * 95) / 4904
  gp1 <- (1009 + 2 * 477 + 3 * 475) / 4904
  hp1 <- hshape / hrate

  n_particles <- nrow(state)

  N <- pars$N
  beta <- pars$beta
  delta <- pars$delta

  ## factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts
  ## per week in long partnerships
  rf <- beta * 1.5 / 7
  rg <- beta * 1 / 7
  mftransm <- MSIf * N * fp1 * rf + MSIg * N * gp1 * rg + beta * MIh * N * hp1

  Y <- rbinom(n_particles, size = ceiling(newI), prob = delta)
  Ytravel <- rbinom(n_particles, size = Y, prob = seedrate / (seedrate + mftransm))
  Yendog <- Y - Ytravel
  Yunk <- 0

  ## c( Ytravel = Ytravel,  Yendog = Y - Ytravel, Yunk = 0)
  ## Y <- Ytravel + Yendog + Yunk - just reconstituting the above

  ## mcstate always works in log space
  log <- TRUE

  ## Erik has this check, but I am not sure it's correct; if we're
  ## working on the entire set of particles it's ok if some are
  ## impossible. But then also not sure how any NA values would get
  ## here, so let's just error instead
  ##
  ## Ah, I see; while rmeas is vectorised, dmeas is not.
  if (any(is.na(Y))) {
    stop("Some NA Y values that need investigation")
  }

  t1 <- rep(-Inf, n_particles)
  i <- newI >= Y
  t1[i] <- dbinom(Y[i], size = ceiling(newI[i]), prob = delta[i], log = TRUE)
  t1[is.na(t1)] <- -Inf

  t2 <- rep(0, n_particles)
  i <- (Ytravel + Yendog) > 0
  t2[i] <- dbinom(Ytravel[i], size = Ytravel[i] + Yendog[i], prob = seedrate / (seedrate + mftransm), log = TRUE)
  t2[is.na(t2)] <- -Inf

  t1 + t2
}


index_m4 <- function(info) {
  list(run = c(newI = info$newI, MSI = info$MSI, MSIg = info$MSIg, MIh = info$MIh, seedrate = info$seedrate),
       state = NULL)
}
