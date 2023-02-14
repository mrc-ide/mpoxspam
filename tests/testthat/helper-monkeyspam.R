evaluate_pgf <- function(f, x) {
  vapply(f, function(fi) fi(x), numeric(1))
}


reference_update_theta_vacc_4.2 <- function(theta_vacc, vacc_amt) {
  f <- function(x) (2809 + 2000 * x + 95 * x^2) / 4904
  g <- function(x) (2943 + 1009 * x + 477 * x^2 + 475 * x^3) / 4904
  hshape <- 0.26
  hrate <- 1.85 * 7
  h <- function(x) (1 - log(x) / hrate)^(-hshape)

  p0 <- f(theta_vacc) * g(theta_vacc) * h(theta_vacc)
  p1 <- max(0.01, p0 - vacc_amt)
  of <- function(lntheta) {
    (f(exp(lntheta)) * g(exp(lntheta)) * h(exp(lntheta)) - p1)^2
  }
  rf <- function(lntheta) {
    f(exp(lntheta)) * g(exp(lntheta)) * h(exp(lntheta)) - p1
  }
  o <- optimise(lower = -1e4, upper = 0, f = of)
  r <- uniroot(lower = -1e4, upper = 0, f = rf)
  c(exp(o$minimum), exp(r$root))
}

reference_update_theta_vacc_4.3 <- function(prop_vacc_targetted) {
  p1 <- max(0.0001, 1 - prop_vacc_targetted)

  hshape <- 0.26
  hrate <- 1.85 * 7
  h <- function(x) (1 - log(x) / hrate)^(-hshape)
  of <- function(lntheta) (h(exp(lntheta)) - p1)^2
  rf <- function(lntheta) h(exp(lntheta)) - p1
  o <- optimise(lower = -1e4, upper = 0, f = of)
  r <- uniroot(lower = -1e4, upper = 0, f = rf)
  c(exp(o$minimum), exp(r$root))
}


reference_compare <- function(state, observed, pars) {
  ## Derived from spam.mpxv@b755562
  log <- TRUE
  delta0 <- pars$delta0
  delta1 <- pars$delta1
  delta_slope <- pars$delta_slope
  Ytravel <- observed$Ytravel
  Yendog <- observed$Yendog
  Yunk <- observed$Yunk

  f <- function(state) {
    time <- state[["time"]]
    newI <- state[["newI"]]
    newIseed <- state[["newIseed"]]

    Y <- ceiling(Ytravel + Yendog + Yunk)
    delta <- max(.01, min(delta1, delta0 + delta_slope * time))
    if (any(is.na(Y))) {
      return(ifelse(log, 0, 1))
    }
    t1 <- ifelse(log, -Inf, 0)
    if (newI >= Y) {
      t1 <- suppressWarnings(dbinom(Y, size = ceiling(newI), prob = delta, log = log))
    } else {
      return(ifelse(log, -Inf, 0))
    }
    if (is.na(t1)) {
      t1 <- ifelse(log, -Inf, 0)
    }

    t2 <- ifelse(log, 0, 1)
    if ((Ytravel + Yendog) > 0) {
      t2 <- dbinom(ceiling(Ytravel), size = ceiling(Ytravel + Yendog), prob = newIseed / newI, log = log)
    }
    if (is.na(t2)) {
      t2 <- ifelse(log, -Inf, 0)
    }

    rv <- ifelse(log, t1 + t2, t1 * t2)
    rv
  }

  apply(state, 2, f)
}


reference_data <- function() {
  dat <- data.frame(
    week = c(11, 14:39),
    Ytravel = c(0, 0, 0, 0, 2, 7, 30, 35, 41, 45, 77, 146, 117, 106, 93, 78,
                109, 86, 54, 40, 33, 19, 20, 12, 16, 5, 14),
    Yendog = c(0, 0, 0, 7, 4, 18, 45, 77, 78, 118, 214, 225, 200, 171, 192,
               143, 143, 140, 108, 34, 42, 36, 27, 17, 25, 10, 13),
    Yunk = 0,
    day = c(14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112,
            119, 126, 133, 140, 147, 154, 161, 168, 175, 182, 189, 196))
  dat$date <- dat$day + as.Date("2022-03-20")
  dat
}


reference_pars <- function() {
  list(beta0 = 1.97848857973857,
       beta_freq = 7,
       beta_sd = 0.15,
       gamma0 = 0.125,
       gamma1 = 0.25,
       etaf = 0.005,
       etag = 0.01,
       N = 750000,
       i0 = 0,
       delta0 = 0.5,
       delta1 = 0.5,
       delta_slope = 0,
       seedrate0 = -0.111546627766835,
       seedrate_sd = 0.75,
       vacc_freq = 1,
       vacc_start_day = 91,
       vacc_fin_day = 146,
       vacc_targetted = 0.8,
       vacc_efficacy = 0.65,
       vacc_doses = 50e3,
       cumulative_partners_days = 90,
       ## New things added by us:
       exp_noise = 1e6,
       kappa_cases = 1,
       rho_travel = 0.5,
       compare_cases = "binom",
       compare_travel = "binom",
       stochastic_behaviour = 1,
       beta_step = c(0, 0),
       dseedrate_step = c(0, 0))
}


filter_data <- function() {
  mcstate::particle_filter_data(reference_data(), "day", 1, 1)
}
