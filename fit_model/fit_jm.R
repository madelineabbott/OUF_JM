################################################################################
# FIT JOINT LONGITUDINAL-SURVIVAL MODEL WITH LATENT ORNSTEIN UHLENBECK PROCESS
################################################################################

library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(survival)

# set working directory
my_wd <- ''

# If using simulated data, specify setting:
setting <- 1 # true parameter values, options: 1 or 2
meas_pattern <- 1 # measurement pattern, options: 1, 2, 3, or 4
g <- 1 # random seed

# Specify density of additional points used to approximate cumul. hazard function
quad_gap <- 0.8 # time between points, manuscript considers 0.2, 0.8, 1.2
# If you don't want additional grid points, use quad_gap = 999

################################################################################
# LOAD DATA
################################################################################

# simulated data:
surv_dat <- read.csv(paste0(my_wd, 'sim_data/surv_dat_s', setting,
                            '_meas_pattern_', meas_pattern, '_g', g, '.csv'))
long_dat <- read.csv(paste0(my_wd, 'sim_data/long_dat_s', setting,
                            '_meas_pattern_', meas_pattern, '_g', g, '.csv'))

N <- nrow(surv_dat) # sample size (number of individuals)

################################################################################
# SET UP GRID POINTS FOR APPROXIMATING THE CUMULATIVE HAZARD FUNCTION
################################################################################

t_max <- max(long_dat$obs_time) # maximum follow-up time
num_quad <- ceiling(t_max / quad_gap) # number of points that we'll add

# set up longer dataframe that has a row for each measurement occasion, each
#  event time, and also each of the added grid points
aug_dat <- data.frame(id = rep(1:N, each = num_quad + 1),
                      time = rep(seq(0, t_max, length.out = num_quad + 1),
                                 times = N),
                      time_type = 'quad')
long_dat_temp <- long_dat %>%
  dplyr::select(id, time) %>%
  mutate(time_type = 'meas')
surv_dat_temp <- surv_dat %>%
  mutate(time = obs_time, time_type = 'event') %>%
  dplyr::select(id, time, time_type)
if (quad_gap == 999){
  # no grid, just measurement occasions and event time
  all_times_template <- rbind(long_dat_temp, surv_dat_temp) %>%
    arrange(id, time, time_type)
}else{
  all_times_template <- rbind(aug_dat, long_dat_temp, surv_dat_temp) %>%
    arrange(id, time, time_type)
}


# remove any grid points after event time
all_times_template <- left_join(all_times_template, surv_dat, by = 'id')
all_times_template <- all_times_template %>%
  filter(time <= obs_time) %>%
  dplyr::select(id, time, time_type)

all_times_template <- all_times_template %>%
  mutate(time = round(time, 10)) %>%
  distinct(id, time, .keep_all = TRUE)

# if the added points are too close to the measurement occasions, then we drop
#  the added points (here we define "too close" as withing 30% of quad_gap)
nrow_orig_template <- nrow(all_times_template)
all_times_template <- all_times_template %>%
  group_by(id) %>%
  mutate(min_gap = abs(pmin(time - lag(time, default = -99),
                            time - lead(time, default = -99)))) %>%
  ungroup() %>%
  mutate(drop = ifelse(time_type == 'quad' & min_gap <= 0.3*quad_gap, 1, 0)) %>%
  filter(drop == 0) %>%
  dplyr::select(-c(min_gap, drop))

################################################################################
# DEFINE VARIABLES FOR STAN
################################################################################

# calculate time between all points (added grid pts + meas. occ. + event times)
all_times_template <- all_times_template %>%
  group_by(id) %>%
  mutate(deltat = time - lag(time, default = 0)) %>%
  ungroup()

# data input for stan:
dat = list(# number of meas. occasions + grid points
           Nall = nrow(all_times_template),
           # number of meas. occasions only (no grid points)
           Nlong = sum(all_times_template$time_type == 'meas'),
           # number of individuals
           N = length(unique(long_dat$id)),
           # number of measured longitudinal outcomes/items
           K = 4,
           # number of latent factors
           P = 2,
           # cumulative # of meas occ. + grid points per person (for generating etas)
           cumu = cumsum(c(table(all_times_template$id))),
           # total # meas occ. + grid points per person (for generating etas)
           repme = c(table(all_times_template$id)), 
           # row numbers corresponding to measurement occasions
           meas_occ_rows = c(1:nrow(all_times_template))[which(
             all_times_template$time_type == 'meas')],
           # cumulative # of meas occ. per person (for observed outcomes Y)
           cumu_meas = cumsum(c(table(all_times_template$id[which(
             all_times_template$time_type == 'meas')]))),
           # total # of meas occ. per person (for observed outcomes Ys)
           repme_meas = c(table(all_times_template$id[which(
             all_times_template$time_type == 'meas')])), 
           # observed longitudinal outcomes Ys (in simulated data, K = 4 outcomes)
           Y = long_dat[,paste0('y', 1:4)],
           # event indicator
           status = surv_dat$status,
           # gap times (for all time points formed by meas. occ. + grid points)
           deltat = all_times_template$deltat)


# initial parameter values
init_params <- list(
  'lambda' = rep(1, 4),
  'sigma_u' = rep(0.1, 4),
  'sigma_eps' = rep(0.1, 4),
  'theta_ou' = array(c(1, 0.5, 0.5, 1), dim = c(2, 2)),
  'rho' = -0.5,
  'beta_0' = 1, 'beta_1' = -1, 'beta_2' = 1)

################################################################################
# FIT MODEL
################################################################################

# fit joint model specified in simulation study (constant baseline hazard)
file <- file.path(paste0(my_wd, 'fit_model/fit_jm_marg_over_u.stan'))

mod <- cmdstan_model(file)

# fit model
fit <- mod$sample(
  data = dat,
  seed = g,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 2000, 
  iter_sampling = 1000,
  init = list(init_params),
  save_warmup = F
)


