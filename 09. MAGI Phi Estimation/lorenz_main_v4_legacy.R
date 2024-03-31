# packages we want to use
library(deSolve)
library(magi)
library(stringr) 
library(glue)

# existing functions that we compressed
source("lorenz_fxns.R")

# let's fix beta
beta = 8/3

#### EXPERIMENTAL SETTINGS (all via command-line) ####

# last Lorenz parameter
rho <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
sigma <- 10.00 # FIXED!

# initial conditions (at t=0.0)
X0 <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
Y0 <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
Z0 <- as.numeric(commandArgs(trailingOnly=TRUE)[4])

# time settings: what interval are we looking at?
t_trunc <- 0.0
t_max <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
  
# mechanical settings. swapped out num_obs for DENSITY of observations!
dens_obs <- as.double(commandArgs(trailingOnly=TRUE)[6])
discretization <- as.double(commandArgs(trailingOnly=TRUE)[7])
alpha <- as.numeric(commandArgs(trailingOnly=TRUE)[8]) # noise level multiplier
specify_noise <- FALSE
hmc_steps <- 101 # irrelevant towards phi estimation, so don't worry!

# final argument: setting a random seed -- the only variable controlled by the main_runscript.sh
seed <- as.double(commandArgs(trailingOnly=TRUE)[9])

# compute the implied num_obs, rounding up to be generous
num_obs <- ceiling( (t_max - t_trunc) * dens_obs )

# need to form our theta vector + initial conditions from inputs
theta <- c(beta, rho, sigma)
x0 <- c(X0, Y0, Z0)

# create a file name to store our logs
fname <- paste("results/", "rho=", rho, "_X0=", X0, "_alpha=", alpha, 
               "_seed=", seed, ".csv", sep="")

#### DATA GENERATION + DISCRETIZATION ####

# generate our true trajectory:NO NOISE! (NON-FORECAST VERSION) + truncate to t_trunc
num_pts <- t_max * 250 + 1
x <- ode(y=x0, times=seq(from=0.0, to=t_max, length.out=num_pts), 
         func=modelODE, parms=theta)
y <- as.data.frame(x)
names(y) <- c("time", "X", "Y", "Z")
y <- y[y$time >= t_trunc,]

# set our seed (also will be our trial number)
set.seed(seed)

# get our observed points that we WILL NOISE + set discretization (EVENLY)
y_obs <- data.frame(y)[seq(from=1, to=nrow(y), by=(nrow(y) %/% num_obs)),]

# get our X0, Y0, Z0 after truncation
X0_true <- y_obs[1,]$X
Y0_true <- y_obs[1,]$Y
Z0_true <- y_obs[1,]$Z

# get our noise levels, for each component: noise = (max - min) * alpha
noise_levels <- as.numeric( (apply(y, 2, max) - apply(y, 2, min) )[-1] * alpha )

# get our noise_draws + combine together
noises_X <- matrix(rnorm(n=nrow(y_obs), mean=0.0, 
                         sd=noise_levels[1]), nrow=nrow(y_obs))
noises_Y <- matrix(rnorm(n=nrow(y_obs), mean=0.0, 
                         sd=noise_levels[2]), nrow=nrow(y_obs))
noises_Z <- matrix(rnorm(n=nrow(y_obs), mean=0.0, 
                         sd=noise_levels[3]), nrow=nrow(y_obs))
noise_draws <- cbind(noises_X, noises_Y, noises_Z)

# add our noise
y_obs[,-1] <- y_obs[,-1] + noise_draws

# set our discretization!
y_obs <- setDiscretization(y_obs, level=discretization)

# generate TRUTH, NO NOISE trajectory, then truncate: "y_true"
num_pts_true <- t_max * 250 + 1
x_true <- ode(y=x0, times=seq(from=0.0, to=t_max, length.out=num_pts_true), 
              func=modelODE, parms=theta)
y_true <- as.data.frame(x_true)
names(y_true) <- c("time", "X", "Y", "Z")
y_true <- y_true[y_true$time >= t_trunc,]

#### DEFINE CONTROLS + RUN FULL STAGE EXPERIMENT ####

# start our list of controls + run standard MAGI, -1 is to omit "time" column.
controls <- list(niterHmc=hmc_steps, burninRatio=0.5)
result <- MagiSolver(y=y_obs[,-1], odeModel=lorenzmodel, 
                     tvec=y_obs$time, control=controls)

#### SAVE OUR RESULTS FOR THE FULL STAGE EXPERIMENT ####
  
# result$phi is where stuff is stored c(result$phi)
new_data <- data.frame(rho=c(rho), X0=c(X0), HMC=c(hmc_steps), TM=c(t_max), 
                       dobs=c(dens_obs), disc=c(discretization), alpha=c(alpha), 
                       SN=c(specify_noise), seed=c(seed), 
                       phi1_X=c(result$phi)[1], phi2_X=c(result$phi)[2], 
                       phi1_Y=c(result$phi)[3], phi2_Y=c(result$phi)[4], 
                       phi1_Z=c(result$phi)[5], phi2_Z=c(result$phi)[6])

# Check if the file exists, create / append new file as necessary
if (file.exists(fname)) {
  write.table(new_data, file = fname, sep = ",", col.names = FALSE, 
              row.names = FALSE, append = TRUE)
} else {
  write.table(new_data, file = fname, sep = ",", col.names = TRUE, 
              row.names = FALSE)
}

