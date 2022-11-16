evaluate_pgf <- function(f, x) {
  vapply(f, function(fi) fi(x), numeric(1))
}


reference_update_theta_vacc <- function(theta_vacc, vacc_amt) {
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
