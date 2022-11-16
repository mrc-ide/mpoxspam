compare_4.2 <- function(state, observed, pars) {
  browser()
  ## Unpack modelled:
  newI <- state["newI", ]
  newIseed <- state["newIseed", ]
  time <- state["time", 1L] # same over all particles

  ## Unpack observed data:
  Ytravel <- observed$Ytravel
  Yendog <- observed$Yendog
  Yunk <- observed$Yunk

  n_particles <- nrow(state)

  ## TODO: need to get time in here, why don't we have this already?
  Y <- ceiling(Ytravel + Yendog + Yunk)
  delta <- max(.01, min(pars$delta1, pars$delta0 + pars$delta_slope * time))
  if (any(is.na(Y))) {
    stop("Unexpected missing value")
  }
  t1 <- rep_len(-Inf, n_particles)
  i <- newI >= Y
  t1[i] <- dbinom(Y, size = ceiling(newI[i]), prob = delta, log = TRUE)
  ## t1[any(is.na(t1))] <- -Inf # hopefully not needed?

  t2 <- rep_len(0, n_particles)
  i <- Ytravel[i] + Yendog[i] > 0
  t2[i] <- dbinom(ceiling(Ytravel[i]),
                  size = ceiling(Ytravel[i] + Yendog[i]),
                  prob = newIseed[i] / (newIseed[i] + newI[i]), log = TRUE)
  ## t2[any(is.na(t2))] <- -Inf # hopefully not needed?

  t1 + t2
}

index_4.2 <- function(info) {
  list(run = c(newI = info$index$new_I,
               newISeed = info$index$newIseed,
               time = info$index$time),
       state = NULL)
}
