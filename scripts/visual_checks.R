library( mpoxspam ) 


reference_pars <- function() {
  list(beta0 = 12,
       beta_freq = 7,
       beta_sd = 0.15,
       gamma0 = 0.125,
       gamma1 = 0.25,
       etaf = 0.005,
       etag = 0.01,
       hshape = 0.26,
       hrate = 12.95,
       N = 750000,
       i0 = 0,
       delta0 = 0.5,
       delta1 = 0.5,
       delta_slope = 0,
       seedrate0 = 0.75, 
       dseedrate0=0,
       seedrate_sd = 0.75,
       vacc_start_day = 91,
       vacc_duration = 55,
       vacc_start_day2 = 91+45,
       vacc_duration2 = 55,
       vacc_targetted = 0.8,
       vacc_efficacy = 0.65,
       vacc_efficacy2 = 1.,
       vacc_doses = 50e3,
       vacc_doses2 = 25e3,
       cumulative_partners_days = 90,
       ## New things added by us:
       exp_noise = 1e6,
       kappa_cases = 1,
       rho_travel = 0.5,
       compare_cases = "binom",
       compare_travel = "binom", 
       dt = .1 ) 
}


pars <- reference_pars()

pars$beta0 <- 6 #12  # high R 4
pars$beta_sd <- 0 # constant transmission rate
pars$N <- 1e5 # small N
# initial burst of seeding only :
pars$dseedrate0 <- -.05
pars$seedrate_sd <- 0
# no vacc 
pars$vacc_start_day = Inf 
pars$vacc_start_day2 = Inf


pars$beta_step <- pars$beta0
pars$dseedrate_step <- pars$dseedrate0
pars$stochastic_behaviour <- 1

set.seed( 1111 )

# remove H contacts
if (F)
{
	pars$hshape = 0.0001
	pars$hrate = 1
}

# remove H contacts & speed up G turnover 
if (F)
{
	pars$hshape = 0.0001
	pars$hrate = 1
	pars$etag = 1/7
}

# heavy vacc on day 50 for one week, all targetted, 100% VE 
if (F)
{
	pars$vacc_start_day = 50
	pars$vacc_start_day2 = Inf
	pars$vacc_efficacy = 1.
	pars$vacc_duration = 7
	pars$vacc_doses = .2*1e5
	pars$vacc_targetted=1.
}

# heavy vacc on day 50 for one week, all random, 100% VE 
if (F)
{
	pars$vacc_start_day = 50
	pars$vacc_start_day2 = Inf
	pars$vacc_efficacy = 1.
	pars$vacc_duration = 7
	pars$vacc_doses = 0.2*1e5
	pars$vacc_targetted=0.
}

# heavy vacc on day 50 for one week, V1 VE 0%, V2 VE 100%, all random
if (T)
{
	pars$vacc_start_day = 50
	pars$vacc_start_day2 = 51
	pars$vacc_efficacy = 0.0
	pars$vacc_efficacy2 = 1.
	pars$vacc_duration = 7
	pars$vacc_duration2 = 7
	pars$vacc_doses = 0.5*1e5
	pars$vacc_doses2 = 0.5*1e5
	pars$vacc_targetted=0.
}

# heavy vacc on day 1 reaching 100% of population, V1 VE 100%, V2 VE 100%, all random
if (F)
{
	pars$vacc_start_day = 1
	pars$vacc_start_day2 = 3
	pars$vacc_efficacy = 1.0
	pars$vacc_efficacy2 = 1.0
	pars$vacc_duration = 1
	pars$vacc_duration2 = 1
	pars$vacc_doses = pars$N-1
	pars$vacc_doses2 = pars$N-1
	pars$vacc_targetted=1e-3
}

# no vacc at all 
if (F)
{
	pars$vacc_start_day = Inf
	pars$vacc_start_day2 = Inf
	pars$vacc_efficacy = 0.0
	pars$vacc_efficacy2 = 0.0
	pars$vacc_duration = 1
	pars$vacc_duration2 = 1
	pars$vacc_doses = 0
	pars$vacc_doses2 = 0
	pars$vacc_targetted=0.
}



# one particle, random seed
m <- model$new(pars, 1, 1
#~ , seed = 20230301
, n_threads = 1) 

tfin <- 400 / 0.1 
taxis <- seq(1, tfin)
res <- m$simulate(taxis)
rownames(res) <- names(m$info()$index)
# last state: 
Y  <- m$state(); rownames(Y) <-  m$info()$index |> names() ; print(Y)  

par( mfrow = c( 3,4 ))
res['I',,] |> plot()
res['R',,] |> plot()
res['V1',,] |> plot()
res['V2',,] |> plot()
#~ res['cutf',,] |> plot()
#~ res['cutg',,] |> plot()
#~ res['cuth',,] |> plot()
#~ res['cuts',,] |> plot()
res['cumulative_partners',,] |> plot()
res['Reff',,] |> plot(log='y'); abline ( h = 1 )
res['Reff_f',,] |> plot(log='y'); abline ( h = 1 )
res['Reff_g',,] |> plot(log='y'); abline ( h = 1 )
res['Reff_h',,] |> plot(log='y'); abline ( h = 1 )
(res['Reff_h',,] + res['Reff_g',,] + res['Reff_f',,]) |> plot(log='y'); abline ( h = 1 )

#~ res['seedrate',,] |> plot()
#~ res['dseedrate',,] |> plot()

#~ X11()
#~ par( mfrow = c(3, 2 ))
#~ res['MSEf',,] |> plot() 
#~ res['MSEg',,] |> plot() 
#~ res['MEh',,] |> plot() 
#~ res['MSIf',,] |> plot() 
#~ res['MSIg',,] |> plot() 
#~ res['MIh',,] |> plot() 

print( 
calcR0(beta = pars$beta0, gamma1 = 1/4, hs = 0.26, hr = 12.95 )
)
