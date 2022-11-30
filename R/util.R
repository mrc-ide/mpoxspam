`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


ll_nbinom <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }

  mu <- model + rexp(length(model), rate = exp_noise)

  dnbinom(data, kappa, mu = mu, log = TRUE)
}


ll_binom <- function(data_x, data_size, model, exp_noise) {
  if (is.na(data_size) || is.na(data_x) ) {
    return(numeric(length(model)))
  }
  prob <- model + rexp(length(model), rate = exp_noise)
  dbinom(data_x, data_size, prob, log = TRUE)
}
