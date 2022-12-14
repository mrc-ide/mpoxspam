##' @importFrom stats rexp
ll_betabinom <- function(data_a, data_b, model_a, model_b, rho, exp_noise) {
  n <- length(model_a)
  if (is.na(data_a) || is.na(data_b) || (data_a == 0 && data_b == 0)) {
    return(numeric(n))
  }
  da <- rexp(n, exp_noise)
  db <- rexp(n, exp_noise)
  ldbetabinom(data_a, data_b, model_a + da, model_b + db, rho)
}

##' @importFrom stats dnbinom
ll_nbinom <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }

  mu <- model + rexp(length(model), rate = exp_noise)

  dnbinom(data, kappa, mu = mu, log = TRUE)
}


ldbetabinom <- function(data_a, data_b, model_a, model_b, rho) {
  prob_a <- model_a / (model_a + model_b)
  prob_b <- model_b / (model_a + model_b)

  a <- prob_a * (1 / rho - 1)
  b <- prob_b * (1 / rho - 1)

  lchoose(data_a + data_b, data_a) +
    lbeta(data_a + a, data_b + b) - lbeta(a, b)
}
