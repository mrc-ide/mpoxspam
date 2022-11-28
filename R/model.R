model_compare <- function(state, observed, pars) {
  ## Unpack modelled:
  newI <- state["newI", ]
  newIseed <- state["newIseed", ]
  time <- state["time", 1L] # same over all particles

  ## Unpack observed data:
  Ytravel <- observed$Ytravel
  Yendog <- observed$Yendog
  Yunk <- observed$Yunk

  n_particles <- ncol(state)

  Y <- ceiling(Ytravel + Yendog + Yunk)
  delta <- max(.01, min(pars$delta1, pars$delta0 + pars$delta_slope * time))
  if (is.na(Y)) {
    return(rep_len(0, n_particles))
  }

  t1 <- rep_len(-Inf, n_particles)
  i <- newI >= Y
  t1[i] <- dbinom(Y, size = ceiling(newI[i]), prob = delta, log = TRUE)

  if (Ytravel + Yendog > 0) {
    i <- newIseed + newI > 0
    t2 <- rep_len(-Inf, n_particles)
    t2[i] <- dbinom(ceiling(Ytravel),
                    size = ceiling(Ytravel + Yendog),
                    prob = newIseed[i] / (newIseed[i] + newI[i]), log = TRUE)
  } else {
    t2 <- 0
  }

  t1 + t2
}

model_index <- function(info) {
  list(run = c(newI = info$index$newI,
               newIseed = info$index$newIseed,
               time = info$index$time),
       state = integer(0))
}


model_filter <- function(data, ..., use_compiled_compare = FALSE) {
  compare <- if (use_compiled_compare) NULL else model_compare
  mcstate::particle_filter$new(data, model, ...,
                               compare = compare,
                               index = model_index)
}