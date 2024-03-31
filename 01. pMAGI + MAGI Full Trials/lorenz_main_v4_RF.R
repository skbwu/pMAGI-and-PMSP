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
sigma <- as.numeric(commandArgs(trailingOnly=TRUE)[2])

# initial conditions (at t=0.0)
X0 <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
Y0 <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
Z0 <- as.numeric(commandArgs(trailingOnly=TRUE)[5])

# time settings: what interval are we looking at?
t_trunc <- as.numeric(commandArgs(trailingOnly=TRUE)[6])
t_max <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
  
# mechanical settings. swapped out num_obs for DENSITY of observations!
dens_obs <- as.double(commandArgs(trailingOnly=TRUE)[8])
discretization <- as.double(commandArgs(trailingOnly=TRUE)[9])
alpha <- as.numeric(commandArgs(trailingOnly=TRUE)[10]) # noise level multiplier
specify_noise <- as.logical(commandArgs(trailingOnly=TRUE)[11])
hmc_steps <- as.double(commandArgs(trailingOnly=TRUE)[12])

# as of 11/24/2023: pilot settings governing length of pilot + no. HMC.
pilot_length <- as.numeric(commandArgs(trailingOnly=TRUE)[13]) # length of our pilot interval
pilot_hmc <- as.double(commandArgs(trailingOnly=TRUE)[14]) # no. of HMC steps for pilot.
pilot_discretization <- as.double(commandArgs(trailingOnly=TRUE)[15]) # no. of HMC steps for pilot.

# final argument: setting a random seed -- the only variable controlled by the main_runscript.sh
seed <- as.double(commandArgs(trailingOnly=TRUE)[16])

# compute the implied num_obs, rounding up to be generous
num_obs <- ceiling( (t_max - t_trunc) * dens_obs )

# need to form our theta vector + initial conditions from inputs
theta <- c(beta, rho, sigma)
x0 <- c(X0, Y0, Z0)

# instantiate + make directory for our pathname
pathname <- paste("results/", "rho=", rho, "_sigma=", sigma,  
                  "_X0=", X0, "_Y0=", Y0, "_Z0=", Z0,
                  "_TT=", t_trunc, "_TM=", t_max, 
                  "_DOBS=", dens_obs, "_D=", discretization,
                  "_alpha=", alpha,  "_SN=", specify_noise, "_HMC=", hmc_steps, 
                  "_PL=", pilot_length, "_PHMC=", pilot_hmc, "_PD=", pilot_discretization,
                  "/seed=", seed, sep="")
dir.create(pathname, recursive=TRUE)

#### DATA GENERATION + DISCRETIZATION ####

# what's the name of the data set we want to use as y_obs?
data_fname_obs <- paste("data/rho=", rho, "_sigma=", sigma, "_X0=", X0, "_Y0=", 
                    Y0, "_Z0=", Z0, "_TT=", t_trunc, "_TM=", t_max, "_DOBS=", 
                    dens_obs, "_alpha=", alpha, "_seed=", seed, 
                    "/y_samp_noised.csv", sep="")
y_obs <- read.csv(data_fname_obs, col.names = c("time", "X", "Y", "Z"))
y_obs <- y_obs[y_obs$time >= t_trunc,]

# set our discretization!
y_obs <- setDiscretization(y_obs, level=discretization)

# what's the name of the data set we want to use for our y_true?
data_fname_interval <- paste("data/rho=", rho, "_sigma=", sigma, "_X0=", X0, "_Y0=", 
                             Y0, "_Z0=", Z0, "_TT=", t_trunc, "_TM=", t_max, "_DOBS=", 
                             dens_obs, "_alpha=", alpha, "_seed=", seed, 
                             "/y_interval.csv", sep="")
y_true <- read.csv(data_fname_interval, col.names = c("time", "X", "Y", "Z"))
y_true <- y_true[y_true$time >= t_trunc,]

#### DEFINE CONTROLS + RUN PILOT STAGE EXPERIMENT to learn phi ####

# check if we're even doing the pilotPhi, based on pilot_length
if (pilot_length > 0.0) {
  
  # make an indicator random variable that we're using pilot
  usePilot <- TRUE
  
  # isolate the data that we are using for pilot: 
  t_pilot_end <- min(t_trunc + pilot_length, t_max)
  y_obs_pilot <- subset(y_obs, time >= t_trunc & time <= t_pilot_end)
  
  # clear the discretizations, then set a custom discretization later.
  y_obs_pilot <- setDiscretization(y_obs_pilot[complete.cases(y_obs_pilot$X),], 
                                   level=pilot_discretization)
  
  # start our list of pilot_controls, respecting useFixedSigma from full trial
  # note that HMC actually matter. We are only getting samples for debugging.
  # i.e., was the phi learned on this smaller interval even a good idea?
  pilot_controls <- list(niterHmc=pilot_hmc, burninRatio=0.5)
  
  # specify noise for the pilot if full-scale requires noise specified, too.
  if (specify_noise == TRUE) {
    pilot_controls$useFixedSigma <- TRUE
    pilot_controls$sigma <- noise_levels
  }
  
  # run our pilot MAGI. The -1 is to omit the "time" column!
  pilot_result <- MagiSolver(y=y_obs_pilot[,-1], odeModel=lorenzmodel, 
                             tvec=y_obs_pilot$time, control=pilot_controls)
  
  # extract out the phi.pilot estimated values
  phi.pilot <- pilot_result$phi
  
  #### SAVE OUR RESULTS FOR DEBUGGING FOR THE PILOT STAGE EXPERIMENT
  
  # save our traces results!
  pilot_traces <- as.data.frame(cbind(pilot_result$theta, pilot_result$lp))
  names(pilot_traces) <- c("beta", "rho", "sigma", "lp")
  write.csv(x=pilot_traces, file=paste(pathname, "pilot_traces.csv", sep="/"), row.names=FALSE)
  
  # save our inferred trajectories + bounds as .csv. start with lower bound.
  pilot_xLB <- apply(pilot_result$xsampled, c(2,3), function(x) quantile(x, 0.025))
  pilot_xLB <- as.data.frame(cbind(y_obs_pilot$time, pilot_xLB))
  names(pilot_xLB) <- c("time", "X", "Y", "Z")
  write.csv(x=pilot_xLB, file=paste(pathname, "pilot_inferred_trajs_LB.csv", sep="/"), row.names=FALSE)
  
  # repeat for the mean
  pilot_xMean <- apply(pilot_result$xsampled, c(2,3), mean)
  pilot_xMean <- as.data.frame(cbind(y_obs_pilot$time, pilot_xMean))
  names(pilot_xMean) <- c("time", "X", "Y", "Z")
  write.csv(x=pilot_xMean, file=paste(pathname, "pilot_inferred_trajs_Mean.csv", sep="/"), row.names=FALSE)
  
  # repeat for the upper bound
  pilot_xUB <- apply(pilot_result$xsampled, c(2,3), function(x) quantile(x, 0.975))
  pilot_xUB <- as.data.frame(cbind(y_obs_pilot$time, pilot_xUB))
  names(pilot_xUB) <- c("time", "X", "Y", "Z")
  write.csv(x=pilot_xUB, file=paste(pathname, "pilot_inferred_trajs_UB.csv", sep="/"), row.names=FALSE)
  
  # save the phi matrices 2x3 dimensions
  filename <- paste(pathname, "pilot_phi.csv", sep="/")
  extract_phi_csv(pilot_result, filename)
  
  # save our inferred parameters as .csv
  filename <- paste(pathname, "pilot_inferred_params.csv", sep="/")
  extract_estimates_csv(pilot_result, filename)
  
  # save our inferred phis as .csv
  filename <- paste(pathname, "pilot_inferred_phis.csv", sep="/")
  extract_estimates_csv(pilot_result, filename)
  
  # SAVE OUR SIGMAS (NOISE)
  pilot_sigmaRecords <- pilot_result$sigma
  write.csv(pilot_sigmaRecords, paste0(pathname, "/pilot_sigmas.csv"))

# else just indicator variable that we're not using the pilot tool.
} else {
  
  # make an indicator random variable that we're NOT using pilot
  usePilot <- FALSE
}


#### DEFINE CONTROLS + RUN FULL STAGE EXPERIMENT ####

# start our list of controls
controls <- list(niterHmc=hmc_steps, burninRatio=0.5)

# make sure to add in our piloted phi into the list of controls, if use pilot.
if (usePilot == TRUE) {
  controls$phi <- phi.pilot
}

# if specifyNoise is True, then specify our TRUE noise levels
if (specify_noise == TRUE) {
  controls$useFixedSigma <- TRUE
  controls$sigma <- noise_levels
  
# fixPhi requires specifying noise. set to pilot-estimated MAP levels.
} else if (usePilot == TRUE) {
    controls$useFixedSigma <- TRUE
    controls$sigma <- apply(pilot_result$sigma, 2, mean)
}

# run our MAGI. The -1 is to omit the "time" column!
result <- MagiSolver(y=y_obs[,-1], odeModel=lorenzmodel, 
                     tvec=y_obs$time, control=controls)

#### SAVE OUR RESULTS FOR THE FULL STAGE EXPERIMENT ####

# save our traces results!
traces <- as.data.frame(cbind(result$theta, result$lp))
names(traces) <- c("beta", "rho", "sigma", "lp")
write.csv(x=traces, file=paste(pathname, "traces.csv", sep="/"), row.names=FALSE)

# save our inferred trajectories + bounds as .csv. start with lower bound.
xLB <- apply(result$xsampled, c(2,3), function(x) quantile(x, 0.025))
xLB <- as.data.frame(cbind(y_obs$time, xLB))
names(xLB) <- c("time", "X", "Y", "Z")
write.csv(x=xLB, file=paste(pathname, "inferred_trajs_LB.csv", sep="/"), row.names=FALSE)
  
# repeat for the mean
xMean <- apply(result$xsampled, c(2,3), mean)
xMean <- as.data.frame(cbind(y_obs$time, xMean))
names(xMean) <- c("time", "X", "Y", "Z")
write.csv(x=xMean, file=paste(pathname, "inferred_trajs_Mean.csv", sep="/"), row.names=FALSE)

# repeat for the upper bound
xUB <- apply(result$xsampled, c(2,3), function(x) quantile(x, 0.975))
xUB <- as.data.frame(cbind(y_obs$time, xUB))
names(xUB) <- c("time", "X", "Y", "Z")
write.csv(x=xUB, file=paste(pathname, "inferred_trajs_UB.csv", sep="/"), row.names=FALSE)
  
# save the phi matrices 2x3 dimensions
filename <- paste(pathname, "phi.csv", sep="/")
extract_phi_csv(result, filename)

# save our inferred parameters as .csv
filename <- paste(pathname, "inferred_params.csv", sep="/")
extract_estimates_csv(result, filename)

# save our inferred phis as .csv
filename <- paste(pathname, "inferred_phis.csv", sep="/")
extract_estimates_csv(result, filename)

# SAVE OUR INITIAL CONDITIONS
y0s <- result$xsampled[,1,]
write.csv(y0s, paste0(pathname, "/init_cond.csv")) 

# SAVE OUR SIGMAS (NOISE)
sigmaRecords <- result$sigma
write.csv(sigmaRecords, paste0(pathname, "/sigmas.csv")) 

# SAVE OUR TRUE PARAMETER VALUES: beta, rho, sigma, noise, initial_conditions: -7.287, -13.919, 11.894
true_params <- t(as.data.frame(c(theta, alpha, t_trunc, t_max)))
colnames(true_params) <- c("beta", "rho", "sigma", "alpha", "t_trunc", "t_max")
write.csv(true_params, paste0(pathname, "/true_params.csv"), row.names=FALSE)

# save our noised observed data, unnoised data (non-forecast), forecasted ground_truth
write.csv(y_obs, paste0(pathname, "/y_obs.csv"), row.names = FALSE)
write.csv(y_true, paste0(pathname, "/y_true.csv"), row.names = FALSE)