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
initial(MIh) <- xinit / 2 #* (1+hpp(1)/hp(1)) / hp(1)
initial(E) <- xinit * N / 2
initial(I) <- xinit * N / 2
initial(newI) <- 0
initial(Eseed) <- 0
initial(newIseed) <- 0
initial(cutf) <- 0
initial(cutg) <- 0
initial(cuth) <- 0
initial(cuts) <- 0
initial(seedrate) <- seedrate0
initial(dseedrate) <- 0
initial(theta_vacc) <- 1
initial(S_vacc) <- 1
initial(beta) <- beta0
initial(cumulative_partners) <- 0
xinit <- i0 / N

## We'll need this; strictly this is (step + 1) * dt but we use unit
## timesteps here.
initial(time) <- step + 1
update(time) <- step + 1

beta0 <- user(2.25)
beta_freq <- user(7)
beta_sd <- user(0.15)  # 1 = 2*sqrt( (75/7)*sigma2 )
gamma0 <- user(0.125)
gamma1 <- user(0.25)
etaf <- user(0.005) ## Anderson Epidemiology 2021
etag <- user(0.01) #Anderson Epidemiology 2021
N <- user(750e3)
## N <- user(N) #
i0 <- user(0)
seedrate0 <- user(0.75)
seedrate_sd <- user(0.75) #2.99 #sd of random walk of daily diff in seedrate
vacc_freq <- user(1)
vacc_start_day <- user(91)
vacc_targetted <- user(.8) # prop vacc targetted vs random
cumulative_partners_days <- user(90)
vacc_duration <- user(55) ## 2022-08-30 - 2022-07-06

## Compare function parameters
delta0 <- user(.50) # ignore.unused
delta1 <- user(.50) # ignore.unused
delta_slope <- user(0.0) # ignore.unused
exp_noise <- user(1e6) # ignore.unused
kappa_cases <- user(5) # ignore.unused
use_nbinom <- user(0, min = 0, max = 1, integer = TRUE) # ignore.unused


## These two were parameters but with derived defaults; could make
## these user again and move that logic into R
vacc_amt <- 0.65 * 50e3 / vacc_duration # assuming about 50k doses by end of august 65% one dose vacc eff
vacc_fin_day <- vacc_start_day + vacc_duration

add_vaccine <-
  (time >= vacc_start_day) &&
  (time <= vacc_fin_day)   &&
  (((time - vacc_start_day) %% vacc_freq) == 0)

amt_targetted <- vacc_targetted * vacc_amt / N
amt_random <- (1 - vacc_targetted) * vacc_amt / N

## Erik's theta_vacc0 is our theta_vacc, his theta_vacc our theta_vacc_use
theta_vacc_use <- if (add_vaccine) update_theta_vacc4_2(theta_vacc, amt_targetted) else theta_vacc
red_f <- (theta_vacc_use * fp(theta_vacc_use)) / (theta_vacc * fp(theta_vacc))
red_g <- (theta_vacc_use * gp(theta_vacc_use)) / (theta_vacc * gp(theta_vacc))

S_vacc_use <- if (add_vaccine) S_vacc * (1 - amt_random) else S_vacc

vaccine_scale_f <- if (add_vaccine) (1 - amt_random) * red_f else 1
vaccine_scale_g <- if (add_vaccine) (1 - amt_random) * red_g else 1

MSEf_vacc <- MSEf * vaccine_scale_f
MSSf_vacc <- MSSf * vaccine_scale_f^2
MSIf_vacc <- MSIf * vaccine_scale_f
MSEg_vacc <- MSEg * vaccine_scale_g
MSSg_vacc <- MSSg * vaccine_scale_g^2
MSIg_vacc <- MSIg * vaccine_scale_g

dot_thetaf <- thetaf * theta_vacc_use
dot_thetag <- thetag * theta_vacc_use
dot_thetah <- thetah * theta_vacc_use

dseedrate_next <- if (time %% beta_freq == 0) rnorm(dseedrate, seedrate_sd) else dseedrate
seedrate_next <- max(as.numeric(0), seedrate + dseedrate_next)
beta_next <- if (time %% beta_freq == 0) max(as.numeric(0), rnorm(beta, beta_sd)) else beta

MSf <- dot_thetaf * S_vacc_use * fp(dot_thetaf) / fp(1)
MSg <- dot_thetag * S_vacc_use * gp(dot_thetag) / gp(1)
# MSh <- dot_thetah * S_vacc_use * hp(dot_thetah) / hp(1) # never used

## factor 1.5 / 7 is act rate per day; factor 1.5 b/c more contacts
## per week in long partnerships
rf <- beta_next * 1.5 / 7
rg <- beta_next * 1 / 7

## TODO: In the pomp code, there was some concern that this would
## become poorly defined and go NA or non-finite
tratef <- max(as.numeric(0), MSIf_vacc * N * fp(1) * rf)
transmf <- rpois(tratef)
dthetaf <- -dot_thetaf * transmf / (MSf * N * fp(1))
dSf <- fp(dot_thetaf) * dthetaf # note prop to transm
meanfield_delta_si_f <- (dot_thetaf * fpp(dot_thetaf) / fp(dot_thetaf))
u1f <- meanfield_delta_si_f
u2f <- (dot_thetaf * fpp(dot_thetaf) + dot_thetaf^2 * fppp(dot_thetaf)) / fp(dot_thetaf)
vf <- u2f - u1f^2
delta_si_f <- (if (transmf == 0) 0
               else rnorm(meanfield_delta_si_f, sqrt(vf / transmf)))

trateg <- max(as.numeric(0), MSIg_vacc * N * gp(1) * rg)
transmg <- rpois(trateg)
dthetag <- -dot_thetag * transmg / (MSg * N * gp(1))
dSg <- gp(dot_thetag) * dthetag # note prop to transm
meanfield_delta_si_g <- (dot_thetag * gpp(dot_thetag) / gp(dot_thetag))
u1g <- meanfield_delta_si_g
u2g <- (dot_thetag * gpp(dot_thetag) + dot_thetag^2 * gppp(dot_thetag)) / gp(dot_thetag)
vg <- u2g - u1g^2
delta_si_g <- (if (transmg == 0) 0
               else rnorm(meanfield_delta_si_g, sqrt(vg / transmg)))

transmseed <- rpois(seedrate_next)

trateh <- max(as.numeric(0), beta_next * MIh * N * hp(1) * S_vacc_use)
transmh <- rpois(trateh) # note S_vacc here, because there is no MSI in MFSH model

dthetah <- -dot_thetah * (transmh + transmseed) / (N * hp(1)) # seeding happens here
dSh <- hp(dot_thetah) * dthetah # note prop to transm
# may also need this separated into seed and non-seed components:
# dSh0 <- hp(dot_thetah)*(-dot_thetah * transmh   / (N*hp(1) ))
# dSh_seed <- hp(dot_thetah)*(-dot_thetah * transmseed  / (N*hp(1) ))
meanfield_delta_si_h <- (1 + dot_thetah * hpp(dot_thetah) / hp(dot_thetah)) # note + 1 for mfsh (vs the above)
u1h <- meanfield_delta_si_h
u2h <- hppp(dot_thetah) * dot_thetah^2 / hp(dot_thetah) + 2 * dot_thetah * hpp(dot_thetah) / hp(dot_thetah) + u1h
vh <- u2h - u1h^2
delta_si_h <- (if (transmh == 0) 0
               else rnorm(meanfield_delta_si_h, sqrt(vh / transmh)))

# record mean and variance of cumulative partners over past x days among new infections
tauf <- (transmf / (transmf + transmg + transmh + transmseed))
taug <- (transmg / (transmf + transmg + transmh + transmseed))
tauh <- ((transmh + transmseed) / (transmf + transmg + transmh + transmseed))
cumulative_partners_next <-
  (1 + etaf * cumulative_partners_days) * (
  tauf * meanfield_delta_si_f +
  taug * dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) +
  tauh * dot_thetaf * fp(dot_thetaf) / f(dot_thetaf)
  ) +
  (1 + etag * cumulative_partners_days) * (
    tauf * dot_thetag * gp(dot_thetag) / g(dot_thetag) +
    taug * meanfield_delta_si_g +
    tauh * dot_thetag * gp(dot_thetag) / g(dot_thetag)
  ) +
  cumulative_partners_days * (
    tauf * dot_thetah * hp(dot_thetah) / h(dot_thetah) +
    taug * dot_thetah * hp(dot_thetah) / h(dot_thetah) +
    tauh * meanfield_delta_si_h
  )

dMSEf <- -gamma1 * MSEf +
  2 * etaf * MSf * MEf -
  etaf * MSEf +
  (-dSf) * (delta_si_f / fp(1)) * (MSSf_vacc / MSf) +
  (-dSg) * (dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) / fp(1)) * (MSSf_vacc / MSf) +
  (-dSh) * (dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) / fp(1)) * (MSSf_vacc / MSf)
dMSIf <- -rf * MSIf -
  gamma1 * MSIf +
  gamma0 * MSEf +
  2 * etaf * MSf * MIf -
  etaf * MSIf +
  (-dSf) * (delta_si_f / fp(1)) * (-MSIf_vacc / MSf) +
  (-dSg) * (dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) / fp(1)) * (-MSIf_vacc / MSf) +
  (-dSh) * (dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) / fp(1)) * (-MSIf_vacc / MSf)

dMSEg <- -gamma0 * MSEg +
  2 * etag * MSg * MEg -
  etag * MSEg +
  (-dSg) * (delta_si_g / gp(1)) * (MSSg_vacc / MSg) +
  (-dSf) * (dot_thetag * gp(dot_thetag) / g(dot_thetag) / gp(1)) * (MSSg_vacc / MSg) +
  (-dSh) * (dot_thetag * gp(dot_thetag) / g(dot_thetag) / gp(1)) * (MSSg_vacc / MSg)
dMSIg <- -rg * MSIg -
  gamma1 * MSIg +
  gamma0 * MSEg +
  2 * etag * MSg * MIg -
  etag * MSIg +
  (-dSg) * (delta_si_g / gp(1)) * (-MSIg_vacc / MSg) +
  (-dSf) * (dot_thetag * gp(dot_thetag) / g(dot_thetag) / gp(1)) * (-MSIg_vacc / MSg) +
  (-dSh) * (dot_thetag * gp(dot_thetag) / g(dot_thetag) / gp(1)) * (-MSIg_vacc / MSg)

dMSSf <- 1 * etaf * MSf^2 - # 2 or 1?
                          etaf * MSSf_vacc -
                          (-dSf) * (delta_si_f / fp(1)) * MSSf_vacc / MSf - # 2 or 1 ?
                          ((-dSg) + (-dSh)) * (dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) / fp(1)) * MSSf_vacc / MSf
dMSSg <- 1 * etag * MSg^2 - # 2 or 1?
                          etag * MSSg_vacc -
                          (-dSg) * (delta_si_g / gp(1)) * MSSg_vacc / MSg - # 2 or 1 ?
                          ((-dSf) + (-dSh)) * (dot_thetag * gp(dot_thetag) / g(dot_thetag) / gp(1)) * MSSg_vacc / MSg

dMEf <- -gamma0 * MEf +
  (-dSf) * (delta_si_f / fp(1)) +
  ((-dSg) + (-dSh)) * (dot_thetaf * fp(dot_thetaf) / f(dot_thetaf) / fp(1))
dMEg <- -gamma0 * MEg +
  (-dSg) * (delta_si_g / gp(1)) +
  ((-dSf) + (-dSh)) * (dot_thetag * gp(dot_thetag) / g(dot_thetag) / gp(1))
dMEh <- -gamma0 * MEh +
  (-dSh) * (delta_si_h / hp(1)) +
  ((-dSf) + (-dSg)) * (dot_thetah * hp(dot_thetah) / h(dot_thetah) / hp(1))

dMIf <- -gamma1 * MIf + gamma0 * MEf
dMIg <- -gamma1 * MIg + gamma0 * MEg
dMIh <- -gamma1 * MIh + gamma0 * MEh

reset_weekly <- step %% 7 == 0

## infected, infectious and not detected:
newE <- transmf + transmg + transmh + transmseed
## exposed, infectious, diagnosed or undiagnosed
I_next <- max(as.numeric(0), I + gamma0 * E - gamma1 * I)
## new case detections. some subset of these will be detected each week
newI_next <- (if (reset_weekly) 0 else newI) + gamma0 * E # will accumulate between obvs
E_next <- max(as.numeric(0), E + newE - gamma0 * E)

# seed state variables
# exposed, infectious, diagnosed or undiagnosed
# new case detections. some subset of these will be detected each week
newIseed_next <- (if (reset_weekly) 0 else newIseed) + gamma0 * Eseed # will accumulate between obvs
Eseed_next <- max(as.numeric(0), Eseed + transmseed - gamma0 * Eseed)

## At this point, Erik checks that there are no na values in
## dtheta[fgh], dMS[EI][fgh], dM[EI][fgh], hopefully we do not need to
## do this.

update(thetaf) <- max(1e-9, min(as.numeric(1), thetaf + dthetaf))
update(MSEf) <- max(as.numeric(0), MSEf_vacc + dMSEf)
update(MEf) <- max(as.numeric(0), MEf + dMEf)
update(MSSf) <- max(as.numeric(0), MSSf_vacc + dMSSf)
update(MSIf) <- max(as.numeric(0), MSIf_vacc + dMSIf)
update(MIf) <- max(as.numeric(0), MIf + dMIf)
update(thetag) <- max(1e-9, min(as.numeric(1), thetag + dthetag))
update(MSEg) <- max(as.numeric(0), MSEg_vacc + dMSEg)
update(MEg) <- max(as.numeric(0), MEg + dMEg)
update(MSSg) <- max(as.numeric(0), MSSg_vacc + dMSSg)
update(MSIg) <- max(as.numeric(0), MSIg_vacc + dMSIg)
update(MIg) <- max(as.numeric(0), MIg + dMIg)
update(thetah) <- max(1e-9, min(as.numeric(1), thetah + dthetah))
update(MEh) <- max(as.numeric(0), MEh + dMEh)
update(MIh) <- max(as.numeric(0), MIh + dMIh)
update(E) <- E_next
update(I) <- I_next
update(newI) <- newI_next
update(Eseed) <- Eseed_next
update(newIseed) <- newIseed_next
update(cutf) <- cutf + transmf
update(cutg) <- cutg + transmg
update(cuth) <- cuth + transmh
update(cuts) <- cuts + transmseed
update(seedrate) <- seedrate_next
update(dseedrate) <- dseedrate_next
update(theta_vacc) <- theta_vacc_use
update(S_vacc) <- S_vacc_use
update(beta) <- beta_next
update(cumulative_partners) <- cumulative_partners_next

config(include) <- "support.hpp"
config(compare) <- "compare.hpp"
