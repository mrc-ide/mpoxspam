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


ll_binom <- function(data_a, data_b, model_a, model_b, exp_noise) {
  if (is.na(data_a) || is.na(data_b) ) {
    return(numeric(length(model)))
  }

  noise_a <- rexp(length(model_a), exp_noise)
  noise_b <- rexp(length(model_b), exp_noise)
  prob <- (model_a + noise_a) / (model_a + model_b + noise_a + noise_b)
  dbinom(data_a, data_a + data_b, prob, log = TRUE)
}


logit <- function(x) {
  log(x / (1 - x))
}
