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
initial( Reff_f ) <- 0 
initial( Reff_g ) <- 0 
initial( Reff_h ) <- 0 
initial( Reff ) <- 0
initial(vacc_period) <- 1
initial(vacc_period2) <- 1

xinit <- i0 / N

## Constants we use in a few places; the as.numeric does a conversion
## to the appropriate floating point type.
fp1 <- fp(as.numeric(1))
gp1 <- gp(as.numeric(1))
hp1 <- hp(as.numeric(1), hshape, hrate)


stochastic_behaviour <- user(1) # Logical switch for user-input trends in beta and seedrate
# only used if stochastic_behaviour == 0, vectors specifying beta and seedrate on each day
beta_step[] <- user()
dseedrate_step[] <- user()
dim(beta_step) <- user()
dim(dseedrate_step) <- user()
dt <- user(.1)
steps_per_week <- 7 / dt

initial(time) <- step
update(time) <- (step + 1) * dt


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
#~ vacc_freq <- user(1)

vacc_targetted <- user(.8) # prop vacc targetted vs random
vacc_efficacy <- user(0.78) # Bertran et el, one dose
vacc_efficacy2 <- user(1.00) # 2 doses, TODO values
cumulative_partners_days <- user(90)

# time-varying doses
vacc_start_day[] <- user()
vacc_start_day2[] <- user()
vacc_doses[] <- user()
vacc_doses2[] <- user() 
vacc_duration[] <- user()
vacc_duration2[] <- user() 

dim(vacc_start_day) <- user()
dim(vacc_start_day2) <- user()
dim(vacc_doses) <- user()
dim(vacc_doses2) <- user() 
dim(vacc_duration) <- user()
dim(vacc_duration2) <- user() 

exp_noise <- user(1e6) # ignore.unused
kappa_cases <- user(1) # ignore.unused
rho_travel <- user(0.5) # ignore.unused
use_new_compare <- user(0) # ignore.unused



# new vacc if within schedule (after delay, taking effect)

vacc_period_inc <- vacc_period + if (time >= vacc_fin_day && (time - dt) < vacc_fin_day) 1 else 0
vacc_period2_inc <- vacc_period2 + if (time >= vacc_fin_day2 && (time - dt) < vacc_fin_day2) 1 else 0

# ensure index does not stray off end of vector into unallocated memory
vacc_period_next <- min(as.integer(vacc_period_inc), length(vacc_start_day))
vacc_period2_next <- min(as.integer(vacc_period2_inc), length(vacc_start_day2))

update(vacc_period) <- vacc_period_next
update(vacc_period2) <- vacc_period2_next
vacc_start_day_step <- vacc_start_day[as.integer(vacc_period)]
vacc_duration_step <- vacc_duration[as.integer(vacc_period)]
vacc_doses_step <- vacc_doses[as.integer(vacc_period)]
vacc_start_day2_step <- vacc_start_day2[as.integer(vacc_period2)]
vacc_duration2_step <- vacc_duration2[as.integer(vacc_period2)]
vacc_doses2_step <- vacc_doses2[as.integer(vacc_period2)]

vacc_fin_day <- vacc_start_day_step + vacc_duration_step
vacc_amt <-  min(vacc_doses_step, N) / (vacc_duration_step/dt)
add_vaccine <- (time >= vacc_start_day_step) && (time < vacc_fin_day)   

vacc_fin_day2 <- vacc_start_day2_step + vacc_duration2_step 
vacc_amt2 <-  min(vacc_doses2_step, N) / (vacc_duration2_step/dt)
add_vaccine2 <- (time >= vacc_start_day2_step) && (time < vacc_fin_day2)   


V1_next <- if (add_vaccine) V1 + vacc_amt / N else V1 
V2_next <- if (add_vaccine2) V2 + vacc_amt2 / N else V2

# effective proportion of population protected by vacc
veff <- if (V1>0) V1*( (V2/V1)*(1-vacc_efficacy)*vacc_efficacy2 + (1-V2/V1)*vacc_efficacy )  else 0
veff_targetted <- veff * vacc_targetted 
veff_untargetted <- veff * (1-vacc_targetted)

theta_vacc_use <- if (add_vaccine) update_theta_vacc4_3(veff_targetted, hshape, hrate) else theta_vacc

vacc_rescale <- if(add_vaccine) 1-(vacc_efficacy*vacc_amt/N + (1-vacc_efficacy)*vacc_efficacy2*vacc_amt2/N) else 1

MSEf_ <- MSEf * vacc_rescale
MSSf_ <- MSSf * vacc_rescale^2
MSIf_ <- MSIf * vacc_rescale
MSEg_ <- MSEg * vacc_rescale
MSSg_ <- MSSg * vacc_rescale^2
MSIg_ <- MSIg * vacc_rescale

MSf <- thetaf * (1-veff) * fp(thetaf) / fp1
MSg <- thetag * (1-veff) * gp(thetag) / gp1
MSh <- (1-veff_untargetted) * thetah*theta_vacc_use* hp(thetah*theta_vacc_use, hshape, hrate) / hp1


## used if stochastic_behaviour == 0
beta_det <- if (as.integer(step) >= length(beta_step))
  beta_step[length(beta_step)] else beta_step[step + 1]
dseedrate_det <- if (as.integer(step) >= length(dseedrate_step))
  dseedrate_step[length(dseedrate_step)] else dseedrate_step[step + 1]

## used if stochastic_behaviour == 1
dseedrate_rw <- if (time %% beta_freq == 0) rnorm(dseedrate, seedrate_sd) else dseedrate
beta_rw <- if (time %% beta_freq == 0) max(as.numeric(0), rnorm(beta, beta_sd)) else beta

dseedrate_next <- if (stochastic_behaviour) dseedrate_rw else dseedrate_det
seedrate_next <- max(as.numeric(0), seedrate + dseedrate_next*dt)
beta_next <- if (stochastic_behaviour) beta_rw else beta_det

## factor 1.5 / 7 is act rate per day; factor 1.5 b/c more contacts
## per week in long partnerships
rf <- beta_next * 1.5 / 7
rg <- beta_next * 1 / 7



tratef <- max(as.numeric(0), ( MSIf_ * N * fp1 * rf) ) 
transmf <- rpois(tratef * dt)
dthetaf <- -thetaf * transmf / (MSf * N * fp1)
#~ dSf <- fp(thetaf) * dthetaf # note prop to transm
dSf <- -transmf / N
meanfield_delta_si_f <- (thetaf * fpp(thetaf) / fp(thetaf))
u1f <- meanfield_delta_si_f
u2f <- (thetaf * fpp(thetaf) + thetaf^2 * fppp(thetaf)) / fp(thetaf)
vf <- min(  u2f - u1f^2, 2*meanfield_delta_si_f )
delta_si_f <- (if (transmf == 0) 0
               else min(meanfield_delta_si_f*2,max(as.numeric(0),rnorm(meanfield_delta_si_f, sqrt(vf / transmf)))))

trateg <- max(as.numeric(0), MSIg_ * N * gp1 * rg)
transmg <- rpois(trateg * dt)
dthetag <- -thetag * transmg / (MSg * N * gp1)
#~ dSg <- gp(thetag) * dthetag # note prop to transm
dSg <- -transmg / N
meanfield_delta_si_g <- (thetag * gpp(thetag) / gp(thetag))
u1g <- meanfield_delta_si_g
u2g <- (thetag * gpp(thetag) + thetag^2 * gppp(thetag)) / gp(thetag)
vg <-  min( u2g - u1g^2, 3*meanfield_delta_si_g)
delta_si_g <- (if (transmg == 0) 0
               else min(meanfield_delta_si_g*2,max(as.numeric(0), rnorm(meanfield_delta_si_g, sqrt(vg / transmg))) ))


dot_thetah <- thetah * theta_vacc_use  

# imported infections
# note imports scaled down by susceptibility in h contacts 
seed_scale_down <- veff_targetted * (dot_thetah * hp(dot_thetah, hshape, hrate) / hp1) + 
  veff_untargetted^2*(thetah * hp(thetah, hshape, hrate) / hp1) + 
  (1-veff)*(thetah * hp(thetah, hshape, hrate) / hp1)
transmseed <- rpois(seedrate_next * dt) * seed_scale_down 

trateh <- beta_next * N * MIh * MSh * hp1
transmh <- rpois(trateh * dt) 


dthetah <- -thetah * (transmh + transmseed) / (N * hp1 * MSh ) # seeding happens here
#~ dSh <- hup(thetah,  vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2 , theta_vacc_use, hshape, hrate) * dthetah *(1-veff) # note prop to transm, (1-veff) added because hu(x) generates normalised dist among unvacc
dSh <- -(transmh + transmseed)/N
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

dMSEf <- -gamma1 * MSEf_ *dt +
   etaf * MSf * MEf*dt - #*2
  etaf * MSEf_*dt +
  (-dSf) * (delta_si_f / fp1) * (MSSf_ / MSf - MSEf_/MSf) +
  (-dSg) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (MSSf_ / MSf - MSEf_/MSf) +
  (-dSh) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (MSSf_ / MSf - MSEf_/MSf)
dMSIf <- -rf * MSIf_*dt -
  gamma1 * MSIf_*dt +
  gamma0 * MSEf_*dt +
   etaf * MSf * MIf*dt - #*2
  etaf * MSIf_*dt +
  (-dSf) * (delta_si_f / fp1) * (-MSIf_ / MSf) +
  (-dSg) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (-MSIf_ / MSf) +
  (-dSh) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * (-MSIf_ / MSf)


dMSEg <- -gamma0 * MSEg_*dt +
   etag * MSg * MEg *dt- #*2 
  etag * MSEg_*dt +
  (-dSg) * (delta_si_g / gp1) * (MSSg_ / MSg) +
  (-dSf) * (thetag * gp(thetag) / g(thetag) / gp1) * (MSSg_ / MSg - MSEg_/MSg) +
  (-dSh) * (thetag * gp(thetag) / g(thetag) / gp1) * (MSSg_ / MSg - MSEg_/MSg)
dMSIg <- -rg * MSIg_*dt -
  gamma1 * MSIg_*dt +
  gamma0 * MSEg_*dt +
   etag * MSg * MIg*dt - #*2
  etag * MSIg_ *dt+
  (-dSg) * (delta_si_g / gp1) * (-MSIg_ / MSg) +
  (-dSf) * (thetag * gp(thetag) / g(thetag) / gp1) * (-MSIg_ / MSg) +
  (-dSh) * (thetag * gp(thetag) / g(thetag) / gp1) * (-MSIg_ / MSg)

dMSSf <- 1 * etaf * MSf^2 *dt- # 2 or 1?
                          etaf * MSSf_*dt -
                          2*(-dSf) * (delta_si_f / fp1) * MSSf_ / MSf - # 2 or 1 ?
                          2*((-dSg) + (-dSh)) * (thetaf * fp(thetaf) / f(thetaf) / fp1) * MSSf_ / MSf
dMSSg <- 1 * etag * MSg^2*dt - # 2 or 1?
                          etag * MSSg_*dt -
                          2*(-dSg) * (delta_si_g / gp1) * MSSg_ / MSg - # 2 or 1 ?
                          2*((-dSf) + (-dSh)) * (thetag * gp(thetag) / g(thetag) / gp1) * MSSg_ / MSg

dMEf <- -gamma0 * MEf *dt +
  (-dSf) * ( (1+delta_si_f) / fp1) +
  ((-dSg) + (-dSh)) * (thetaf * fp(thetaf) / f(thetaf) / fp1)
dMEg <- -gamma0 * MEg*dt +
  (-dSg) * ( (1+delta_si_g) / gp1) +
  ((-dSf) + (-dSh)) * (thetag * gp(thetag) / g(thetag) / gp1)


hup1 <- hup(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) 
hu1 <- hu(thetah,vacc_targetted, V1, V2, vacc_efficacy, vacc_efficacy2, theta_vacc_use, hshape, hrate) 
dMEh <- -gamma0 * MEh*dt +
  (-dSh) * (delta_si_h / hp1) + # +1 for delta_si_h already included 
  ((-dSf) + (-dSg)) * (thetah * hup1 / hu1 / hp1 )


dMIf <- -gamma1 * MIf*dt + gamma0 * MEf*dt
dMIg <- -gamma1 * MIg*dt + gamma0 * MEg*dt
dMIh <- -gamma1 * MIh*dt + gamma0 * MEh*dt

reset_weekly <- step %% steps_per_week == 0

## infected, infectious and not detected:
newE <- transmf + transmg + transmh + transmseed
## exposed, infectious, diagnosed or undiagnosed
I_next <- max(as.numeric(0), I + gamma0 * E*dt - gamma1 * I*dt)
## new case detections. some subset of these will be detected each week
newI_next <- (if (reset_weekly) 0 else newI) + gamma0 * E*dt # will accumulate between obvs
E_next <- max(as.numeric(0), E + newE - gamma0 * E*dt)
S_next <- max(as.numeric(0), S - newE)
R_next <- max(as.numeric(0), R + gamma1 * I*dt)

# seed state variables
# exposed, infectious, diagnosed or undiagnosed
# new case detections. some subset of these will be detected each week
newIseed_next <- (if (reset_weekly) 0 else newIseed) + gamma0 * Eseed *dt# will accumulate between obvs
Eseed_next <- max(as.numeric(0), Eseed + transmseed - gamma0 * Eseed*dt)


## Need as.numeric here to get a good cast to the correct real type
update(thetaf) <- max(as.numeric(1e-9), min(as.numeric(1), thetaf + dthetaf))
update(MSEf) <- max(as.numeric(0), min(as.numeric(1), MSEf_ + dMSEf) )
update(MEf) <- max(as.numeric(0), MEf + dMEf)
update(MSSf) <- max(as.numeric(0), min( as.numeric(1), MSSf_ + dMSSf) )
update(MSIf) <- max(as.numeric(0), min( as.numeric(1) ,  MSIf_ + dMSIf) )
update(MIf) <- max(as.numeric(0), MIf + dMIf)
update(thetag) <- max(as.numeric(1e-9), min(as.numeric(1), thetag + dthetag))
update(MSEg) <- max(as.numeric(0), min( as.numeric(1), MSEg_ + dMSEg))
update(MEg) <- max(as.numeric(0), MEg + dMEg)
update(MSSg) <- max(as.numeric(0), min( as.numeric(1), MSSg_ + dMSSg) )
update(MSIg) <- max(as.numeric(0), min( as.numeric(1), MSIg_ + dMSIg) )
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

update(Reff_f) <- tratef / I / gamma1
update(Reff_g) <- trateg  / I / gamma1
update(Reff_h) <- trateh / I / gamma1
update( Reff ) <- (tratef + trateg + trateh) / I / gamma1


config(include) <- "support.hpp"
config(compare) <- "compare.hpp"

# debugging 
#~ print( "time: {time; .0f} N: {N} MIh {MIh}  MSh {MSh} hp1 {hp1} trateh {trateh}", when= trateh > 30)
#~ print('transmf {transmf} transmg {transmg} transmh {transmh} transmseed {transmseed}', when= trateh > 30)
#~ print( 'dSh {dSh} dSf {dSf} dSg {dSg} delta_si_h {delta_si_h} MEh {MEh}' , when= trateh > 30)
#~ print( 'thetah {thetah} hup1 {hup1} hu1 {hu1} ' , when= trateh > 30)
#~ print( 'MSEf {MSEf}  dMSEf {dMSEf} MSEf_ {MSEf_}' , when = trateh > 30 )
#~ print( 'MSIf {MSIf}  dMSIf {dMSIf} MSIf_ {MSIf_}' , when = trateh > 30 )
#~ print("time: {time; .0f} veff: {veff} trateh {trateh} dt are you there? {dt}", when= trateh > 30)
#~ print( 'MSh {MSh} ' )
#~ print( '{fp1} {gp1} {hp1} ' )
#~ print ( '..........................{1}', when = trateh > 30 )
# print( "time: {time} V1: {V1} V2: {V2} vacc_amt {vacc_amt} vacc_amt2 {vacc_amt2} vacc_fin_day {vacc_fin_day} vacc_fin_day2 {vacc_fin_day2} vacc_period: {vacc_period} vacc_duration_step {vacc_duration_step}")
# print( "time: {time} step: {tmp} length: {tmp2} dseedrate: {dseedrate_det} seedrate {seedrate} seedrate_next {seedrate_next} dseedrate_next {dseedrate_next}")


# print("veff: {veff} tratef: {tratef} trateg {trateg} trateh: {trateh} seedrate_next {seedrate_next} V1: {V1_next} V2: {V2_next} R: {R} beta_next {beta_next} MSf: {MSf} dMSIf: {dMSIf}")



