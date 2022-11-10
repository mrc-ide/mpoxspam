initial(thetaf) <- 1
initial(MSEf) <- 0
initial(MEf) <- 0
initial(MSSf) <- 1
initial(MSIf) <- 0
initial(MIf) <- 0
initial(thetag) <- 1
initial(MSEg) <- 0
initial(MEg) <- 0
initial(MSSg) <- 1
initial(MSIg) <- 0
initial(MIg) <- 0
initial(thetah) <- 1 - xinit
initial(MEh) <- xinit / 2
initial(MIh) <- xinit / 2 # (1 + hpp(1) / hp1) / hp1
initial(E) <- xinit * N / 2
initial(I) <- xinit * N / 2
initial(newI) <- 0
initial(cutf) <- 0
initial(cutg) <- 0
initial(cuth) <- 0
initial(cuts) <- 0
initial(seedrate) <- seedrate0
initial(theta_vacc) <- 1

xinit <- i0 / N

## Modified on entry
##
## E

## Modified on exit
##
## thetaf
## MEf
## MIf
## thetag
## MEg
## MIg
## thetah
## MEh
## MIh

## Modified both
##
## MSEf
## MSSf
## MSIf
## MSEg
## MSSg
## MSIg

time <- step
add_vaccine <-
  (time >= vacc_start_day) &&
  (time <= vacc_fin_day)   &&
  (((time - vacc_start_day) %% vacc_freq) == 0)

vaccine_scale <- if (add_vaccine) 1 - vacc_amt else 1
MSEf_vacc <- MSEf * vaccine_scale
MSSf_vacc <- MSSf * vaccine_scale^2
MSIf_vacc <- MSIf * vaccine_scale
MSEg_vacc <- MSEg * vaccine_scale
MSSg_vacc <- MSSg * vaccine_scale^2
MSIg_vacc <- MSIg * vaccine_scale

dot_thetah <- thetah * theta_vacc # was .thetah
seedrate_next <- max(0, seedrate + rnorm(0, seedrate_sd))

## Fraction of contacts in set S in the f, g, and h partnership networks
MSf <- thetaf * pgf_f1 / fp1
MSg <- thetag * pgf_g1 / gp1
## This quantity is unused, including in the original version
# MSh <- dot_thetah * pgf_h1 / hp1

## beta is transmission probability per contact * contact rate
## inflation factor. Transmission rates:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5142082/
##
## factor 1.5 / 7 is sex act rate per day; factor 1.5 b/c more
## contacts per week in long partnerships
rf <- beta * 1.5 / 7
rg <- beta * 1 / 7

transmf <- rpois(MSIf_vacc * N * fp1 * rf)
dthetaf <- -thetaf * transmf / (MSf * N * fp1)
dSf <- pgf_f1 * dthetaf # note prop to transm
meanfield_delta_si_f <- (thetaf * pgf_f2 / pgf_f1)
u1f <- meanfield_delta_si_f
u2f <- (thetaf * pgf_f2 + thetaf^2 * pgf_f3) / pgf_f1
vf <- u2f - u1f^2
delta_si_f <- (if (transmf == 0) 0
               else rnorm(meanfield_delta_si_f, sqrt(vf / transmf)))

transmg <- rpois(MSIg_vacc * N * gp1 * rg)
dthetag <- -thetag * transmg / (MSg * N * gp1)
dSg <- pgf_g1 * dthetag # note prop to transm
meanfield_delta_si_g <- (thetag * pgf_g2 / pgf_g1)
u1g <- meanfield_delta_si_g
u2g <- (thetag * pgf_g2 + thetag^2 * pgf_g3) / pgf_g1
vg <- u2g - u1g^2
delta_si_g <- (if (transmg == 0) 0
               else rnorm(meanfield_delta_si_g, sqrt(vg / transmg)))

transmseed <- rpois(seedrate_next)
transmh <- rpois(beta * MIh * N * hp1)
dthetah <- -dot_thetah * (transmh + transmseed) / (N * hp1) # seeding happens here
dSh <- pgf_h1 * dthetah # note prop to transm
meanfield_delta_si_h <- (1 + dot_thetah * pgf_h2 / pgf_h1) # note + 1 for mfsh
u1h <- meanfield_delta_si_h
u2h <- pgf_h3 * dot_thetah^2 / pgf_h1 + 2 * dot_thetah * pgf_h2 / pgf_h1 + u1h
vh <- u2h - u1h^2
delta_si_h <- (if (transmh == 0) 0
               else rnorm(meanfield_delta_si_h, sqrt(vh / transmh)))

## TODO:
##
## Remove double +'s
## Find common factors, especially in last two lines of each
## Structure with arrays: (fgh) first, then SI then SI
dMSEf <- -gamma1 * MSEf_vacc +
  +2 * etaf * MSf * MEf +
   -etaf * MSEf_vacc +
    +(-dSf) * (delta_si_f / fp1) * (MSSf_vacc / MSf) +
    +(-dSg) * (thetaf * pgf_f1 / pgf_f0 / fp1) * (MSSf_vacc / MSf) +
    +(-dSh) * (thetaf * pgf_f1 / pgf_f0 / fp1) * (MSSf_vacc / MSf)
dMSIf <- -rf * MSIf_vacc +
  -gamma1 * MSIf_vacc +
   +gamma0 * MSEf_vacc +
    +2 * etaf * MSf * MIf +
     -etaf * MSIf_vacc +
      +(-dSf) * (delta_si_f / fp1) * (-MSIf_vacc / MSf) +
      +(-dSg) * (thetaf * pgf_f1 / pgf_f0 / fp1) * (-MSIf_vacc / MSf) +
      +(-dSh) * (thetaf * pgf_f1 / pgf_f0 / fp1) * (-MSIf_vacc / MSf)

dMSEg <- -gamma0 * MSEg_vacc +
  +2 * etag * MSg * MEg +
   -etag * MSEg_vacc +
    +(-dSg) * (delta_si_g / gp1) * (MSSg_vacc / MSg) +
    +(-dSf) * (thetag * pgf_g1 / pgf_g0 / gp1) * (MSSg_vacc / MSg) +
    +(-dSh) * (thetag * pgf_g1 / pgf_g0 / gp1) * (MSSg_vacc / MSg)
dMSIg <- -rg * MSIg_vacc +
  -gamma1 * MSIg_vacc +
   +gamma0 * MSEg_vacc +
    +2 * etag * MSg * MIg +
     -etag * MSIg_vacc +
      +(-dSg) * (delta_si_g / gp1) * (-MSIg_vacc / MSg) +
      +(-dSf) * (thetag * pgf_g1 / pgf_g0 / gp1) * (-MSIg_vacc / MSg) +
      +(-dSh) * (thetag * pgf_g1 / pgf_g0 / gp1) * (-MSIg_vacc / MSg)

dMSSf <- +1 * etaf * MSf^2 +
  -etaf * MSSf_vacc +
  -(-dSf) * (delta_si_f / fp1) * MSSf_vacc / MSf +
  -((-dSg) + (-dSh)) * (thetaf * pgf_f1 / pgf_f0 / fp1) * MSSf_vacc / MSf
dMSSg <- +1 * etag * MSg^2 +
  -etag * MSSg_vacc +
  -(-dSg) * (delta_si_g / gp1) * MSSg_vacc / MSg +
  -((-dSf) + (-dSh)) * (thetag * pgf_g1 / pgf_g0 / gp1) * MSSg_vacc / MSg

dMEf <- -gamma0 * MEf +
  +(-dSf) * (delta_si_f / fp1) +
  +((-dSg) + (-dSh)) * (thetaf * pgf_f1 / pgf_f0 / fp1)
dMEg <- -gamma0 * MEg +
  +(-dSg) * (delta_si_g / gp1) +
  +((-dSf) + (-dSh)) * (thetag * pgf_g1 / pgf_g0 / gp1)
dMEh <- -gamma0 * MEh +
  +(-dSh) * (delta_si_h / hp1) +
  +((-dSf) + (-dSg)) * (dot_thetah * pgf_h1 / pgf_h0 / hp1)

dMIf <- -gamma1 * MIf + gamma0 * MEf
dMIg <- -gamma1 * MIg + gamma0 * MEg
dMIh <- -gamma1 * MIh + gamma0 * MEh

## infected, infectious and not detected:
newE <- transmf + transmg + transmh + transmseed
## exposed, infectious, diagnosed or undiagnosed
E_next <- max(0, E + newE - gamma0 * E)
I_next <- max(0, I + gamma0 * E - gamma1 * I)

update(thetaf) <- max(1e-9, min(1, thetaf + dthetaf))
update(MSEf) <- max(0, MSEf_vacc + dMSEf)
update(MEf) <- max(0, MEf + dMEf)
update(MSSf) <- max(0, MSSf_vacc + dMSSf)
update(MSIf) <- max(0, MSIf_vacc + dMSIf)
update(MIf) <- max(0, MIf + dMIf)
update(thetag) <- max(1e-9, min(1, thetag + dthetag))
update(MSEg) <- max(0, MSEg_vacc + dMSEg)
update(MEg) <- max(0, MEg + dMEg)
update(MSSg) <- max(0, MSSg_vacc + dMSSg)
update(MSIg) <- max(0, MSIg_vacc + dMSIg)
update(MIg) <- max(0, MIg + dMIg)
update(thetah) <- max(1e-9, min(1, thetah + dthetah))
update(MEh) <- max(0, MEh + dMEh)
update(MIh) <- max(0, MIh + dMIh)
update(E) <- E_next
update(I) <- I_next
update(newI) <- newI + gamma0 * E_next
## cumulative transmission
update(cutf) <- cutf + transmf
update(cutg) <- cutg + transmg
update(cuth) <- cuth + transmh
update(cuts) <- cuts + transmseed
update(seedrate) <- seedrate_next
update(theta_vacc) <- theta_vacc

## Parameters
beta <- user(2.25)
gamma0 <- user(1 / 8)
gamma1 <- user(1 / 4)
etaf <- user(1 / 200) ## Anderson Epidemiology 2021
etag <- user(1 / 100) # Anderson Epidemiology 2021
N <- user(.49 * 56.3 * (1 - .24) * 1e6 * 0.02) # male*england(millions)*adult*msm
i0 <- user(0)
## 40pc detection rate
delta <- user(.4) # ignore.unused
seedrate0 <- user(0.75)
seedrate_sd <- user(0.05)
vacc_freq <- user(7)
vacc_amt <- user(0.04)
vacc_start_day <- user(91)
vacc_fin_day <- user(91 + 5 * 7)

fp1 <- (2000 + 2 * 95) / 4904
gp1 <- (1009 + 2 * 477 + 3 * 475) / 4904
hp1 <- hshape / hrate

pgf_f0 <- (2809 + 2000 * thetaf + 95 * thetaf^2) / 4904
pgf_f1 <- (2000 + 2 * 95 * thetaf) / 4904
pgf_f2 <- (2 * 95) / 4904
pgf_f3 <- 0

pgf_g0 <- (2943 + 1009 * thetag + 477 * thetag^2 + 475 * thetag^3) / 4904
pgf_g1 <- (1009 + 2 * 477 * thetag + 3 * 475 * thetag^2) / 4904
pgf_g2 <- (2 * 477 + 2 * 3 * 475 * thetag) / 4904
pgf_g3 <- (2 * 3 * 475) / 4904

hshape <- 0.26
hrate <- 1.85 * 7

pgf_h0 <- (1 - log(dot_thetah) / hrate)^(-hshape)
pgf_h1 <- hshape * (1 - log(dot_thetah) / hrate)^(-hshape - 1) / ((hrate * dot_thetah))
pgf_h2 <- hshape * ((hrate - log(dot_thetah)) / hrate)^(-hshape) * (-hrate + hshape + log(dot_thetah) + 1) / (dot_thetah^2 * (hrate - log(dot_thetah))^2)
pgf_h3 <- hshape * ((hrate - log(dot_thetah)) / hrate)^(-hshape) * (hshape^2 + 3 * hshape + 2 * (hrate - log(dot_thetah))^2 - 3 * (hrate - log(dot_thetah)) * (hshape + 1) + 2) / (dot_thetah^3 * (hrate - log(dot_thetah))^3)
