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
##' @importFrom stats dbinom
model_compare <- function(state, observed, pars) {

  ## Unpack modelled:
  n_particles <- ncol(state)
  time <- state["time", 1L] # same over all particles
  delta <- max(.01, min(pars$delta1, pars$delta0 + pars$delta_slope * time))

  newI <- state["newI", ] ## total new cases (E->I)
  newIseed <- state["newIseed", ] ## new travel-associated cases
  newIendog <- newI - newIseed ## new locally-acquired cases
  ## allowing for under-ascertainment
  model_Y <- ceiling(newI * delta)
  model_Ytravel <- ceiling(newIseed * delta)
  model_Yendog <- model_Y - model_Ytravel

  ## Unpack observed data:
  Ytravel <- observed$Ytravel
  Yendog <- observed$Yendog
  Yunk <- observed$Yunk
  Y <- ceiling(Ytravel + Yendog + Yunk)

  if (is.na(Y)) {
    return(rep_len(0, n_particles))
  }

  if (pars$compare_cases == "negbinom") {
    # Yendog ~ NegBin(model_Yendog, kappa)
    ll_cases <- ll_nbinom(Yendog, model_Yendog, pars$kappa_model_Y,
                          pars$exp_noise)
  } else if(pars$compare_cases == "binom") {
    # Y ~ Bin(model_Y, delta)
    ll_cases <- rep_len(-Inf, n_particles)
    i <- newI >= Y
    ll_cases[i] <- dbinom(Y, size = ceiling(newI[i]), prob = delta, log = TRUE)
  } else {
    stop(sprintf("unrecognised compare function %s", pars$compare_cases))
  }

  if (pars$compare_travel == "betabinom") {
    # Ytravel ~ BetaBinom(Y, model_p, rho)
    ll_travel <- ll_betabinom(Ytravel, Yendog, newIseed, newIendog,
                              pars$rho_travel, pars$exp_noise)
  } else if (pars$compare_travel == "binom") {
    # Ytravel ~ Bin(Y, model_p)
    if (Ytravel + Yendog > 0) {
      i <- newIseed + newI > 0
      ll_travel <- rep_len(-Inf, n_particles)
      ll_travel[i] <- dbinom(ceiling(Ytravel),
                             size = ceiling(Ytravel + Yendog),
                             prob = newIseed[i] / newI[i],
                             log = TRUE)
    } else {
      ll_travel <- 0
    }
  } else if (pars$compare_travel == "negbinom") {
    # Y ~ NegBin(model_Ytravel * delta, kappa)
    ll_travel <- ll_nbinom(Ytravel, model_Ytravel, pars$kappa_travel,
                           pars$exp_noise)
  } else {
    stop(sprintf("unrecognised compare function %s", pars$compare_travel))
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
##' @param generator Alternative generator to use, you'll need this if
##'   using a GPU version.
##'
##' @export
model_filter <- function(data, ..., generator = NULL,
                         use_compiled_compare = FALSE) {
  compare <- if (use_compiled_compare) NULL else model_compare
  mcstate::particle_filter$new(data, model = generator %||% model, ...,
                               compare = compare,
                               index = model_index)
}
