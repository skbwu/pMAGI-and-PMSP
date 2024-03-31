# packages we want to use
library(deSolve)
library(magi)
library(stringr) 
library(glue)

# existing functions that we compressed
source("lorenz_fxns.R")

# let's fix beta
beta = 8/3
sigma <- 10.0

#### EXPERIMENTAL SETTINGS (all via command-line) ####

# trajectory-specifying parameters
rho <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
X0 <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
Y0 <- X0
Z0 <- Y0

# time settings: what interval are we looking at?
TT <- 0.0
TM <- as.numeric(commandArgs(trailingOnly=TRUE)[3]) # {2, 4, 6, 8, 10}
specify_noise <- as.logical(commandArgs(trailingOnly=TRUE)[4])
hmc_pilot_steps <- 4001 # no. of HMC steps for pilot.
  
# mechanical settings
DOBS <- as.double(commandArgs(trailingOnly=TRUE)[5]) # {10, 20, 40}
in_discretization <- as.double(commandArgs(trailingOnly=TRUE)[6]) # {1, 2}
out_discretization <- in_discretization
alpha <- as.numeric(commandArgs(trailingOnly=TRUE)[7]) # noise level multiplier 0.15 or 0.000015

# as of 11/24/2023: pilot settings governing length of pilot + no. HMC.
pilot_length <- 1.0 # length of our pilot interval
pilot_discretization <- 1 # this seems to be the best setting!

# need to form our theta vector + initial conditions from inputs
theta <- c(beta, rho, sigma)
x0 <- c(X0, Y0, Z0)

# forecast-specific settings
TFC <- 5.0 # how far ahead are we forecasting?
forecast_stepsize <- as.numeric(commandArgs(trailingOnly=TRUE)[8])
WSIO <- "last-ode"

# one last setting for how we are doing our pilot work with sequential forecasting
pilot_method <- c("no_pilot", "front_pilot", "back_pilot", "front_online_pilot", "back_online_pilot")[strtoi(commandArgs(trailingOnly=TRUE)[9])]

# final argument: setting a random seed
seed <- as.double(commandArgs(trailingOnly=TRUE)[10])

# to keep apples-to-apples, we will change HMC steps based on stepsize
if (forecast_stepsize == 0.5) {
  hmc_init_steps <- 16001
  hmc_peak_steps <- 16001
} else if (forecast_stepsize == 1.0) {
  hmc_init_steps <- 32001
  hmc_peak_steps <- 32001
} else {
  hmc_init_steps <- 128001
  hmc_peak_steps <- 128001 
}

# instantiate + make directory for our pathname
pathname <- paste("results/", "rho=", rho,
                  "_X0=", X0, "_TM=", TM, 
                  "_DOBS=", DOBS, "_D=", in_discretization,
                  "_alpha=", alpha,  "_SN=", specify_noise, 
                  "_PM=", gsub("_", "-", pilot_method),
                  "_FS=", forecast_stepsize, 
                  "/seed=", seed, sep="")
dir.create(pathname, recursive=TRUE)


#### DATA GENERATION + DISCRETIZATION ####

# a) load in our noisy y_obs, b) cutoff to our obs. range, c) discretize.
data_fname <- paste("forecast_data/rho=", rho, "_sigma=", sigma, "_X0=", X0, "_Y0=", 
                    Y0, "_Z0=", Z0, "_TT=0_TM=10_TFC=5_DOBS=40_alpha=", alpha, "_seed=", seed, 
                    "/y_samp+forecast_noised.csv", sep="")
y_obs_full <- read.csv(data_fname, col.names = c("time", "X", "Y", "Z"))
y_obs <- y_obs_full[y_obs_full$time <= TM,]

# ENFORCE THE DOBS RESTRICTIONS!
if (DOBS != 40) {
  y_obs <- y_obs[((seq(length(y_obs$time)) %% (40 / DOBS)) - 1) == 0,] 
}

# set the discretizations USING THE in-discretization!
y_obs <- setDiscretization(y_obs, level=in_discretization)

# enable forecasting by padding TFC units of time with NA, splice together, RESPECTING THE OUT DISCRETIZATION!
y_add <- subset(y_obs_full, time > TM & time <= TM + TFC) # y_obs_full is always at 40 obs / unit time!
y_add <- setDiscretization(y_add, level=out_discretization)
y_add[, !(names(y_add) %in% c("time"))] <- NaN
y_obs_forecast <- rbind(y_obs, y_add)

# b) load the ground truth too just to see what our true sigmas are
y_true <- read.csv(paste("forecast_data/rho=", rho, "_sigma=", sigma, "_X0=", X0, "_Y0=", 
                         Y0, "_Z0=", Z0, "_TT=0_TM=10_TFC=5_DOBS=40_alpha=", alpha, "_seed=", seed, 
                         "/y_samp+forecast_unnoised.csv", sep=""), col.names = c("time", "X", "Y", "Z"))
sigma_true <- as.numeric((apply(y_true[,-1], 2, max) - apply(y_true[,-1], 2, min)) * alpha)

#### DEFINE CONTROLS + RUN PILOT STAGE EXPERIMENT, IF CALLED FOR to learn phi ####

# are we even doing a pilot? if not, specify that we are not
if (grepl("no_pilot", pilot_method)) {
  
  # indicator variable
  usePilot <- FALSE
  
} else {
  
  # set our indicator variable
  usePilot <- TRUE
  
  # start our list of pilot_controls, respecting useFixedSigma from full trial
  # note that HMC actually matter. We are only getting samples for debugging.
  # i.e., was the phi learned on this smaller interval even a good idea?
  pilot_controls <- list(niterHmc=hmc_pilot_steps, burninRatio=0.5)
  
  # check if we need to add in the true sigma
  if (specify_noise == TRUE) {
    pilot_controls$useFixedSigma <- TRUE
    pilot_controls$sigma <- sigma_true
  }
  
  # are we doing front or back pilot? if front pilot ...
  if (grepl("front", pilot_method)) {
    
    # isolate the data that we are using for FRONT pilot (TT, TT+PL): 
    t_pilot_end <- min(TT + pilot_length, TM)
    y_obs_pilot <- subset(y_obs, time >= TT & time <= t_pilot_end)
  
  # else, we're doing back pilot: (TM-PL, TM)
  } else {
    
    # do our BACK pilot setup
    t_pilot_start <- max(TT, TM - pilot_length)
    y_obs_pilot <- subset(y_obs, time >= t_pilot_start & time <= TM)
    
  }
  
  # clear the discretizations, then set a custom discretization later.
  y_obs_pilot <- setDiscretization(y_obs_pilot[complete.cases(y_obs_pilot$X),], 
                                   level=pilot_discretization)
  
  # run our pilot MAGI. The -1 is to omit the "time" column!
  pilot_result <- MagiSolver(y=y_obs_pilot[,-1], odeModel=lorenzmodel, 
                             tvec=y_obs_pilot$time, control=pilot_controls)
  
  # extract out the phi.pilot estimated values
  phi.pilot <- pilot_result$phi
  
  #### SAVE OUR RESULTS FOR DEBUGGING FOR THE PILOT STAGE EXPERIMENT ####
  
  # save the phi matrices 2x3 dimensions
  filename <- paste(pathname, "pilot_phi.csv", sep="/")
  extract_phi_csv(pilot_result, filename)
  
  # save our inferred parameters as .csv
  filename <- paste(pathname, "pilot_inferred_params.csv", sep="/")
  extract_estimates_csv(pilot_result, filename)
  
  # SAVE OUR SIGMAS (NOISE)
  pilot_sigmaRecords <- pilot_result$sigma
  write.csv(pilot_sigmaRecords, paste0(pathname, "/pilot_sigmas.csv"))
}


#### DEFINE CONTROLS + RUN FULL STAGE FORECASTING ####

# how many forecasting steps do we need?
num_forecast_steps <- TFC / forecast_stepsize

# iterate through each of our forecasting steps
for ( step_num in seq(num_forecast_steps) ) {
  
  # 0. start our timer
  # start_time <- as.numeric(Sys.time())
  
  # 1. for this step, how far are we forecasting into the future?
  t_forecast_seq <- TM + (step_num * forecast_stepsize)
  
  # 1.5 start our list of controls
  controls <- list()
  
  # 2. do we have any sort of pilot module at play?
  if (grepl("no_pilot", pilot_method)) {
    
    # all we need to check is if we need to specify the TRUE noise
    if (specify_noise == TRUE) {
      controls$useFixedSigma <- TRUE
      controls$sigma <- sigma_true
    }
  
  # if yes pilot, do we need to run a new pilot or keep the old pilot?
  } else {
    
    # do we just keep the old pilot or do we have to learn a new pilot?
    if (grepl("online", pilot_method)) {
      
      # if not the first step, we need to learn a new pilot!
      if (step_num == 1) {
        
        # keep our old pilot!
        controls$phi <- phi.pilot
        controls$useFixedSigma <- TRUE
        controls$sigma <- apply(pilot_result$sigma, 2, mean)
        
      } else {
        
        # a) need to run a new pilot! Should have xMean available from prev. step!
        seq_pilot_controls <- list()
        seq_pilot_controls$niterHmc <- 101

        # b) always doing back-pilot 
        t_pilot_end_seq <- xMean$time[length(xMean$time)]
        t_pilot_start_seq <- max(TT, t_pilot_end_seq - forecast_stepsize)
        y_obs_pilot_seq <- subset(xMean, time > t_pilot_start_seq & time <= t_pilot_end_seq)
        
        # c) run our pilot again
        seq_pilot_result <- MagiSolver(y=y_obs_pilot_seq[,-1], odeModel=lorenzmodel, 
                                       tvec=y_obs_pilot_seq$time, control=seq_pilot_controls)
        
        # d) load up our phi parameters into the full controls
        controls$phi <- seq_pilot_result$phi
        controls$useFixedSigma <- TRUE
        controls$sigma <- apply(pilot_result$sigma, 2, mean) # ORIGINAL PILOT!
      }
      
      
    # just keep our original pilot from before
    } else {
      
      # noise should also already be encoded in our pilot
      controls$phi <- phi.pilot
      controls$useFixedSigma <- TRUE
      controls$sigma <- apply(pilot_result$sigma, 2, mean)
      
    }
  }
  
  # 3.1 separate into two cases: 1st forecasting step or not?
  if (step_num == 1) {
    
    # a) go with the specified command-line arguments, more burn-in!
    controls$niterHmc <- hmc_init_steps
    controls$burninRatio <- 0.5
    
    # b) warm-start theta using pilot, if the pilot exists!
    if (!grepl("no_pilot", pilot_method)) {
      controls$thetaInit <- colMeans(pilot_result$theta) 
    }
    
  # 3.2 NOT first forecasting step  
  } else {
    
    # a) compute HMC prop. to peak steps, less burn-in because warm start!
    controls$niterHmc <- ceiling((t_forecast_seq / (TFC + TM)) * hmc_peak_steps)
    controls$burninRatio <- 0.2
    
    # b) set up initialization for theta + setup initialization for xInit_seq
    # note that we are ONLY using "last-ode"! also, mean-ode" is NOT the same meaning as before!
      
    # i) take the last sampled theta as our initialization
    controls$thetaInit <- result_seq$theta[nrow(result_seq$theta),]
    
    # ii) use the last inferred trajs. as init. in-sample
    xInit_seq_prev <- result_seq$xsampled[dim(result_seq$xsampled)[1],,]
    
    # iii) warm start out-of-sample xInit_seq using last(theta), last(y0)
    theta_oos_init <- result_seq$theta[nrow(result_seq$theta),]
    y0_oos_init <- result_seq$xsampled[dim(result_seq$xsampled)[1], 
                                       dim(result_seq$xsampled)[2],]
    
    # c) numerically integrate to get xInit_seq_new
    t_end_prev <- result_seq$tvec[length(result_seq$tvec)]
    ts_oos_init <- y_obs_forecast[y_obs_forecast$time >= t_end_prev 
                                  & y_obs_forecast$time <= t_forecast_seq,]$time
    xInit_seq_new <- ode(y=y0_oos_init, times=ts_oos_init, 
                         func=modelODE, parms=theta_oos_init)[-1,-1]
    
    # d) assemble our concatenated xInit_seq + add to our controls
    xInit_seq <- rbind(xInit_seq_prev, xInit_seq_new)
    controls$xInit <- xInit_seq
  }
  
  # 4. what is our new y_obs input into pMAGI at this step?
  y_obs_seq <- y_obs_forecast[y_obs_forecast$time <= t_forecast_seq,]
  
  # 5. run pMAGI with all inputs readied. omit the time column
  result_seq <- MagiSolver(y=y_obs_seq[,-1], odeModel=lorenzmodel, 
                           tvec=y_obs_seq$time, control=controls)
  
  #### 6. save traces + forecasts for this sequential step, status update ####
  pathname_seq <- paste0(pathname, "/step=", step_num)
  dir.create(pathname_seq, recursive=TRUE)
  
  # save our traces results!
  traces <- as.data.frame(cbind(result_seq$theta, result_seq$lp))
  names(traces) <- c("beta", "rho", "sigma", "lp")
  write.csv(x=traces, file=paste(pathname_seq, "traces.csv", sep="/"), row.names=FALSE)
  
  # save our inferred trajectories + bounds as .csv. start with lower bound.
  xLB <- apply(result_seq$xsampled, c(2,3), function(x) quantile(x, 0.025))
  xLB <- as.data.frame(cbind(y_obs_seq$time, xLB))
  names(xLB) <- c("time", "X", "Y", "Z")
  write.csv(x=xLB, file=paste(pathname_seq, "inferred_trajs_LB.csv", sep="/"), row.names=FALSE)
  
  # repeat for the mean
  xMean <- apply(result_seq$xsampled, c(2,3), mean)
  xMean <- as.data.frame(cbind(y_obs_seq$time, xMean))
  names(xMean) <- c("time", "X", "Y", "Z")
  write.csv(x=xMean, file=paste(pathname_seq, "inferred_trajs_Mean.csv", sep="/"), row.names=FALSE)
  
  # repeat for the upper bound
  xUB <- apply(result_seq$xsampled, c(2,3), function(x) quantile(x, 0.975))
  xUB <- as.data.frame(cbind(y_obs_seq$time, xUB))
  names(xUB) <- c("time", "X", "Y", "Z")
  write.csv(x=xUB, file=paste(pathname_seq, "inferred_trajs_UB.csv", sep="/"), row.names=FALSE)
  
  # -1. end our timer
  # end_time <- as.numeric(Sys.time())
  
  # status update
  # print(paste("Finished step ", step_num, " of ", num_forecast_steps, " in ", end_time - start_time, " seconds.", sep=""))
}
