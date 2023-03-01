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
initial(S) <- N - i0
initial(E) <- i0 / 2
initial(I) <- i0 / 2
initial(R) <- 0
initial(newI) <- 0
initial(Eseed) <- 0
initial(newIseed) <- 0
initial(cutf) <- 0
initial(cutg) <- 0
initial(cuth) <- 0
initial(cuts) <- 0
initial(seedrate) <- seedrate0
initial(dseedrate) <- dseedrate0
initial(theta_vacc) <- 1
initial(beta) <- beta0
initial(cumulative_partners) <- 0
initial( V1 ) <- 0 
initial( V2 ) <- 0


xinit <- i0 / N

## Constants we use in a few places; the as.numeric does a conversion
## to the appropriate floating point type.
fp1 <- fp(as.numeric(1))
gp1 <- gp(as.numeric(1))
hp1 <- hp(as.numeric(1), hshape, hrate)


## We'll need this; strictly this is (step + 1) * dt but we use unit
## timesteps here.
initial(time) <- step + 1
update(time) <- step + 1

stochastic_behaviour <- user(1) # Logical switch for user-input trends in beta and seedrate
# only used if stochastic_behaviour == 0, vectors specifying beta and seedrate on each day
beta_step[] <- user()
dseedrate_step[] <- user()
dim(beta_step) <- user()
dim(dseedrate_step) <- user()

beta0 <- user(2.25)
beta_freq <- user(7)
beta_sd <- user(0.15)  # 1 = 2*sqrt( (75/7)*sigma2 )
gamma0 <- user(0.125)
gamma1 <- user(0.25)
etaf <- user(0.005) ## Anderson Epidemiology 2021
etag <- user(0.01) #Anderson Epidemiology 2021
hshape <- user(0.26) ## Weiss 2020
hrate <- user(12.95) ## Weiss 2020 1.85 * 7 converted to weekly
N <- user(750e3)
i0 <- user(0)
delta0 <- user(.50) # ignore.unused
delta1 <- user(.50) # ignore.unused
delta_slope <- user(0.0) # ignore.unused
seedrate0 <- user(0.75)
dseedrate0 <- user(0)
seedrate_sd <- user(0.75) #2.99 #sd of random walk of daily diff in seedrate
vacc_freq <- user(1)
vacc_start_day <- user(91)
vacc_start_day2 <- user(91+45) #TODO find date 2nd dose provided
vacc_targetted <- user(.8) # prop vacc targetted vs random
vacc_efficacy <- user(0.78) # Bertran et el, one dose
vacc_efficacy2 <- user(1.00) # 2 doses, TODO values
vacc_doses <- user(50e3) # Assume 50k doses in UK
vacc_doses2 <- user(15e3) # TODO find number of 2nd doses
cumulative_partners_days <- user(90)
vacc_duration <- user(55) ## 2022-08-30 - 2022-07-06
vacc_duration2 <- user(55) ## TODO find dates 2nd dose provided

exp_noise <- user(1e6) # ignore.unused
kappa_cases <- user(1) # ignore.unused
rho_travel <- user(0.5) # ignore.unused
use_new_compare <- user(0) # ignore.unused

# new vacc if within schedule (after delay, taking effect)
vacc_amt <-  vacc_doses / (vacc_duration/vacc_freq)
vacc_amt2 <- vacc_doses2 / (vacc_duration2/vacc_freq) # 2nd dose
vacc_fin_day <- vacc_start_day + vacc_duration
vacc_fin_day2 <- vacc_start_day2 + vacc_duration2


add_vaccine <-
  (time >= vacc_start_day) &&
  (time <= vacc_fin_day)   &&
  (((time - vacc_start_day) %% vacc_freq) == 0)

add_vaccine2 <-
  (time >= vacc_start_day2) &&
  (time <= vacc_fin_day2)   &&
  (((time - vacc_start_day2) %% vacc_freq) == 0)



V1_next <- if (add_vaccine) V1 + vacc_amt / N else V1 
V2_next <- if ( add_vaccine2 ) V2 + vacc_amt2 / N else V2
# effective proportion of population protected by vacc
veff <- if (V1>0) V1*( (V2/V1)*(1-vacc_efficacy)*vacc_efficacy2 + (1-V2/V1)*vacc_efficacy )  else 0
#~ S_vacc_use <- if (add_vaccine) S_vacc * (1 - amt_random - amt_targetted) else S_vacc

theta_vacc_use <- if (add_vaccine) update_theta_vacc4_3(veff, hshape, hrate) else theta_vacc

vaccine_scale_f <- if(add_vaccine) 1-(vacc_efficacy*vacc_amt/N + (1-vacc_efficacy)*vacc_efficacy2*vacc_amt2/N) else 1
vaccine_scale_g <- vaccine_scale_f

MSEf_vacc <- MSEf * vaccine_scale_f
MSSf_vacc <- MSSf * vaccine_scale_f^2
MSIf_vacc <- MSIf * vaccine_scale_f
MSEg_vacc <- MSEg * vaccine_scale_g
MSSg_vacc <- MSSg * vaccine_scale_g^2
MSIg_vacc <- MSIg * vaccine_scale_g

MSf <- thetaf * (1-veff) * fp(thetaf) / fp1
MSg <- thetag * (1-veff) * gp(thetag) / gp1

## used if stochastic_behaviour == 0
beta_det <- if (as.integer(step) >= length(beta_step))
  beta_step[length(beta_step)] else beta_step[step + 1]
dseedrate_det <- if (as.integer(step) >= length(dseedrate_step))
  dseedrate_step[length(dseedrate_step)] else dseedrate_step[step + 1]

## used if stochastic_behaviour == 1
dseedrate_rw <- if (time %% beta_freq == 0) rnorm(dseedrate, seedrate_sd) else dseedrate
beta_rw <- if (time %% beta_freq == 0) max(as.numeric(0), rnorm(beta, beta_sd)) else beta

dseedrate_next <- if (stochastic_behaviour) dseedrate_rw else dseedrate_det
seedrate_next <- max(as.numeric(0), seedrate + dseedrate_next)
beta_next <- if (stochastic_behaviour) beta_rw else beta_det

## factor 1.5 / 7 is act rate per day; factor 1.5 b/c more contacts
## per week in long partnerships
rf <- beta_next * 1.5 / 7
rg <- beta_next * 1 / 7



## become poorly defined and go NA or non-finite
tratef <- max(as.numeric(0), ( MSIf_vacc * N * fp1 * rf) )
transmf <- rpois(tratef)
dthetaf <- -thetaf * transmf / (MSf * N * fp1)
dSf <- fp(thetaf) * dthetaf # note prop to transm
meanfield_delta_si_f <- (thetaf * fpp(thetaf) / fp(thetaf))
u1f <- meanfield_delta_si_f
u2f <- (thetaf * fpp(thetaf) + thetaf^2 * fppp(thetaf)) / fp(thetaf)
vf <- min(  u2f - u1f^2, 2*meanfield_delta_si_f )
delta_si_f <- (if (transmf == 0) 0
               else min(meanfield_delta_si_f*2,max(as.numeric(0),rnorm(meanfield_delta_si_f, sqrt(vf / transmf)))))

trateg <- max(as.numeric(0), MSIg_vacc * N * gp1 * rg)
transmg <- rpois(trateg)
dthetag <- -thetag * transmg / (MSg * N * gp1)
dSg <- gp(thetag) * dthetag # note prop to transm
meanfield_delta_si_g <- (thetag * gpp(thetag) / gp(thetag))
u1g <- meanfield_delta_si_g
u2g <- (thetag * gpp(thetag) + thetag^2 * gppp(thetag)) / gp(thetag)
vg <-  min( u2g - u1g^2, 3*meanfield_delta_si_g)
delta_si_g <- (if (transmg == 0) 0
               else min(meanfield_delta_si_g*2,max(as.numeric(0), rnorm(meanfield_delta_si_g, sqrt(vg / transmg))) ))

transmseed <- rpois(seedrate_next) * (thetah * hp(thetah, hshape, hrate) / hp1)


trateh <- max(as.numeric(0), beta_next * MIh * N * hp1 * (S/N)*(1-V1*vacc_efficacy+V2*(1-vacc_efficacy)*vacc_efficacy2) )
transmh <- rpois(trateh) 



dot_thetah <- thetah * theta_vacc_use

MSh <- thetah * hup(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) / hp1 # only unvacc
dthetah <- -thetah * (transmh + transmseed) / (N * hp1 * MSh ) # seeding happens here
dSh <- hup(thetah,  vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2 , theta_vacc_use, hshape, hrate) * dthetah *(1-veff) # note prop to transm, (1-veff) added because hu(x) generates normalised dist among unvacc
meanfield_delta_si_h <- (1 + thetah * hupp(thetah, vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) /
                           hup(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate)) # note + 1 for mfsh (vs f & g)
u1h <- meanfield_delta_si_h
u2h <- huppp(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) * thetah^2 /
  hup(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2,theta_vacc_use, hshape, hrate) +
  2 * thetah * hupp(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2,theta_vacc_use, hshape, hrate) /
  hup(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2,theta_vacc_use, hshape, hrate) + u1h
vh <- min( u2h - u1h^2 , 2*meanfield_delta_si_h)
delta_si_h <- (if (transmh == 0) 0
               else min(meanfield_delta_si_h*2,max(as.numeric(0),rnorm(meanfield_delta_si_h, sqrt(vh / transmh)))))

# record mean and variance of cumulative partners over past x days among new infections
tauf <- (transmf / (transmf + transmg + transmh + transmseed))
taug <- (transmg / (transmf + transmg + transmh + transmseed))
tauh <- ((transmh + transmseed) / (transmf + transmg + transmh + transmseed))
cumulative_partners_next <-
  (1 + etaf * cumulative_partners_days) * (
  tauf * meanfield_delta_si_f +
  taug * thetaf * fp(thetaf) / f(thetaf) +
  tauh * thetaf * fp(thetaf) / f(thetaf)
  ) +
  (1 + etag * cumulative_partners_days) * (
    tauf * thetag * gp(thetag) / g(thetag) +
    taug * meanfield_delta_si_g +
    tauh * thetag * gp(thetag) / g(thetag)
  ) +
  cumulative_partners_days * (
    tauf * dot_thetah * hp(dot_thetah, hshape, hrate) /
      h(dot_thetah, hshape, hrate) +
    taug * dot_thetah * hp(dot_thetah, hshape, hrate) /
      h(dot_thetah, hshape, hrate) +
    tauh * meanfield_delta_si_h
  )

dMSEf <- -gamma1 * MSEf +
  2 * etaf * MSf * MEf -
  etaf * MSEf +
  (-dSf) * (delta_si_f / fp1) * (MSSf_vacc / MSf) +
  (-dSg) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (MSSf_vacc / MSf) +
  (-dSh) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (MSSf_vacc / MSf)
dMSIf <- -rf * MSIf -
  gamma1 * MSIf +
  gamma0 * MSEf +
  2 * etaf * MSf * MIf -
  etaf * MSIf +
  (-dSf) * (delta_si_f / fp1) * (-MSIf_vacc / MSf) +
  (-dSg) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (-MSIf_vacc / MSf) +
  (-dSh) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (-MSIf_vacc / MSf)

dMSEg <- -gamma0 * MSEg +
  2 * etag * MSg * MEg -
  etag * MSEg +
  (-dSg) * (delta_si_g / gp1) * (MSSg_vacc / MSg) +
  (-dSf) * (thetag * gp(thetag) / g(thetag) / gp1) * (MSSg_vacc / MSg) +
  (-dSh) * (thetag * gp(thetag) / g(thetag) / gp1) * (MSSg_vacc / MSg)
dMSIg <- -rg * MSIg -
  gamma1 * MSIg +
  gamma0 * MSEg +
  2 * etag * MSg * MIg -
  etag * MSIg +
  (-dSg) * (delta_si_g / gp1) * (-MSIg_vacc / MSg) +
  (-dSf) * (thetag * gp(thetag) / g(thetag) / gp1) * (-MSIg_vacc / MSg) +
  (-dSh) * (thetag * gp(thetag) / g(thetag) / gp1) * (-MSIg_vacc / MSg)

dMSSf <- 1 * etaf * MSf^2 - # 2 or 1?
                          etaf * MSSf_vacc -
                          (-dSf) * (delta_si_f / fp1) * MSSf_vacc / MSf - # 2 or 1 ?
                          ((-dSg) + (-dSh)) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * MSSf_vacc / MSf
dMSSg <- 1 * etag * MSg^2 - # 2 or 1?
                          etag * MSSg_vacc -
                          (-dSg) * (delta_si_g / gp1) * MSSg_vacc / MSg - # 2 or 1 ?
                          ((-dSf) + (-dSh)) * (thetag * gp(thetag) / g(thetag) / gp1) * MSSg_vacc / MSg

dMEf <- -gamma0 * MEf +
  (-dSf) * (delta_si_f / fp1) +
  ((-dSg) + (-dSh)) * (thetaf * fp(thetaf) / f(thetaf) / fp1)
dMEg <- -gamma0 * MEg +
  (-dSg) * (delta_si_g / gp1) +
  ((-dSf) + (-dSh)) * (thetag * gp(thetag) / g(thetag) / gp1)
dMEh <- -gamma0 * MEh +
  (-dSh) * (delta_si_h / hp1) + 
  ((-dSf) + (-dSg)) * (thetah * hup(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) / hu(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) / hp1 )

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
S_next <- max(as.numeric(0), S - newE)
R_next <- max(as.numeric(0), R + gamma1 * I)

# seed state variables
# exposed, infectious, diagnosed or undiagnosed
# new case detections. some subset of these will be detected each week
newIseed_next <- (if (reset_weekly) 0 else newIseed) + gamma0 * Eseed # will accumulate between obvs
Eseed_next <- max(as.numeric(0), Eseed + transmseed - gamma0 * Eseed)


## Need as.numeric here to get a good cast to the correct real type
update(thetaf) <- max(as.numeric(1e-9), min(as.numeric(1), thetaf + dthetaf))
update(MSEf) <- max(as.numeric(0), min(as.numeric(1), MSEf_vacc + dMSEf) )
update(MEf) <- max(as.numeric(0), MEf + dMEf)
update(MSSf) <- max(as.numeric(0), min( as.numeric(1), MSSf_vacc + dMSSf) )
update(MSIf) <- max(as.numeric(0), min( as.numeric(1) ,  MSIf_vacc + dMSIf) )
update(MIf) <- max(as.numeric(0), MIf + dMIf)
update(thetag) <- max(as.numeric(1e-9), min(as.numeric(1), thetag + dthetag))
update(MSEg) <- max(as.numeric(0), min( as.numeric(1), MSEg_vacc + dMSEg))
update(MEg) <- max(as.numeric(0), MEg + dMEg)
update(MSSg) <- max(as.numeric(0), min( as.numeric(1), MSSg_vacc + dMSSg) )
update(MSIg) <- max(as.numeric(0), min( as.numeric(1), MSIg_vacc + dMSIg) )
update(MIg) <- max(as.numeric(0), MIg + dMIg)
update(thetah) <- max(as.numeric(1e-9), min(as.numeric(1), thetah + dthetah))
update(MEh) <- max(as.numeric(0), MEh + dMEh)
update(MIh) <- max(as.numeric(0), MIh + dMIh)
update(S) <- S_next
update(E) <- E_next
update(I) <- I_next
update(R) <- R_next
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
update(beta) <- beta_next
update(cumulative_partners) <- cumulative_partners_next
update(V1) <- V1_next
update(V2) <- V2_next


config(include) <- "support.hpp"
config(compare) <- "compare.hpp"
