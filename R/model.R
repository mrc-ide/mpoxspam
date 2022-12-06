##' @name model
##' @title Pair approximation transmission model of Mpox virus
##'
##' @description The basic model; we may add more or adapt this one
##'   over time.
##'
##' @export model
NULL

##' Compare observed and modelled data from the basic model. This
##' conforms to the mcstate interface.
##'
##' @title Compare observed and modelled data for the basic model
##'
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [model_index()] so contains
##'   rows corresponding to new infections and scaled time.
##'
##' @param observed Observed data. This will be a list with elements
##'   `Ytravel`, `Yendog` and `Yunk`.
##'
##' @param pars A list of parameters
##'
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##'
##' @export
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

  if (use_nbinom) {
    ll_cases <- ll_nbinom(Y, cases, kappa = pars$kappa_cases, pars$exp_noise)
    ll_travel <- ll_binom(Ytravel, Yendog, newIseed, newI, pars$exp_noise)
  } else {
    ll_cases <- rep_len(-Inf, n_particles)
    i <- newI >= Y
    ll_cases[i] <- dbinom(Y, size = ceiling(newI[i]), prob = delta, log = TRUE)

    if (Ytravel + Yendog > 0) {
      i <- newIseed + newI > 0
      ll_travel <- rep_len(-Inf, n_particles)
      ll_travel[i] <- dbinom(
        ceiling(Ytravel),
        size = ceiling(Ytravel + Yendog),
        prob = newIseed[i] / (newIseed[i] + newI[i]), log = TRUE)
    } else {
      ll_travel <- 0
    }
  }
  ll_cases + ll_travel
}


##' Index of "interesting" elements for the model. This function
##' conforms to the mcstate interface.
##'
##' @title Index of model variables
##'
##' @param info The result of running the `$info()` method on an
##'   initialised [model]
##'
##' @return A list with element `run`, indicating the locations of the
##'   compartments used in [model_compare]
##'
##' @export
model_index <- function(info) {
  list(run = c(newI = info$index$newI,
               newIseed = info$index$newIseed,
               time = info$index$time),
       state = c(cumulative_partners = info$index$cumulative_partners,
                 seedrate = info$index$seedrate,
                 I = info$index$I,
                 newI = info$index$newI,
                 newIseed = info$index$newIseed))
}


##' Create a particle filter; just a convenience function around some
##' mcstate defaults.
##'
##' @title Create particle filter
##'
##' @param data Data to use
##'
##' @param ... Additional arguments to [mcstate::particle_filter]'s constructor
##'
##' @param use_compiled_compare Logical, indicating if we should use
##'   the compiled compare function rather than [model_compare]. This
##'   will likely be faster if several threads are used.
##'
##' @export
model_filter <- function(data, ..., use_compiled_compare = FALSE) {
  compare <- if (use_compiled_compare) NULL else model_compare
  mcstate::particle_filter$new(data, model, ...,
                               compare = compare,
                               index = model_index)
}
