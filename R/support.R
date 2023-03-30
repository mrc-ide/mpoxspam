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


#' @export
pgf_f <- function(x) {
  .Call(`_mpoxspam_test_f`, x)
}

#' @export 
pgf_g <- function(x) {
  .Call(`_mpoxspam_test_g`, x)
}

#' @export
pgf_h <- function(x, hs, hr) {
  .Call(`_mpoxspam_test_h`, x, hs, hr)
}




#' @export 
calcR0 <- function(beta = 10, gamma1 = 1/4, hs = 0.26, hr = 12.95 )
{
	rf <- beta * 1.5/7# factor 1.5/7 is act rate per day; factor 1.5 b/c more contacts per week in long partnerships 
	rg <- beta * 1/7
	pf = rf / (rf + gamma1)
	pg = rg / (rg + gamma1)
	
	hpp <- pgf_h(1, hs, hr)[3]
	hp <- pgf_h(1, hs, hr)[2] 
	gpp <- pgf_g(1)[3]
	gp <- pgf_g(1)[2] 
	fpp <- pgf_f(1)[3]
	fp <- pgf_f(1)[2] 
	
	#NOTE no eta terms b/c of separate timescales
	Rff = pf*fpp/fp 
	Rgg = pg*gpp/gp
	Rhh = (1+hpp/hp) * beta / gamma1
	Rgf = Rhf = pf * fp 
	Rfg = Rhg = pg * gp 
	Rfh = Rgh = (beta / gamma1)*hp 
	K = matrix( c( 
		  Rff, Rfg, Rfh
		, Rgf, Rgg, Rgh
		, Rhf, Rhg, Rhh 
	) , nrow = 3, byrow=TRUE )
	R0 = eigen(K)$values[1] 
	c( R0 = R0 , fR0 = Rff, gR0 = Rgg, hR0 = Rhh )
}
