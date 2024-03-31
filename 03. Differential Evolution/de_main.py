import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math, sys, copy, os, shutil
from scipy.optimize import differential_evolution
from scipy.integrate import solve_ivp

# fixed settings
seeds, TT, alpha = 10, 0.0, 0.15

# what are our possible regime settings? (all bound to sys.argv[1]!)
regime = ["stable_canonical", "stable_transient-chaos", 
          "chaotic_butterfly", "chaotic_no-butterfly"][int(sys.argv[1])]
theta_true = [(8/3, 6.0, 10.0), (8/3, 23.0, 10.0), 
              (8/3, 28.0, 10.0), (8/3, 28.0, 10.0)][int(sys.argv[1])]
y0 = [(5.0, 5.0, 5.0), (5.0, 5.0, 5.0), 
      (5.0, 5.0, 5.0), (2.0, 2.0, 2.0)][int(sys.argv[1])]

# settings to iterate over (separate cmd-line arguments)
dobs = [5, 10, 20, 40][int(sys.argv[2])]
TM = [2.0, 4.0, 6.0, 8.0, 10.0][int(sys.argv[3])]

# Differential Evolution method-specific parameters: CONSULT PROF. YANG ON HOW IN-DEPTH TO GRID-SEARCH!
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html

# multiplier for setting total population size 
# (default is 15, but since we have 3 vars, let's go 45)
popsizes = [45]

# should be between 0-2, can specify as range for random oscillation per generation
mutations = [(0.5, 1.0), (0.5, 1.5), (1.0, 1.5)]

# should be between 0-1, controls mutants vs. population stability.
recombinations = [0.7, 0.5, 0.3]

# create a directory to store our results
if "results" not in os.listdir(): os.mkdir("results")
    
    
# make a directory for this regime
if regime not in os.listdir("results"): os.mkdir(f"results/{regime}")
    
# extract our true beta, rho, sigma
beta_true, rho_true, sigma_true = theta_true

# create the name of the .csv file that we are going to use for this set of trials.
FNAME = f"TM={TM}_alpha={alpha}_dobs={dobs}.csv"

# create our dataframe to store our results for this Lorenz sim, cols = PSO settings.
df = pd.DataFrame(data=None, columns=["popsize", "mut_LB", "mut_UB", "recomb", 
                                      "seed", "beta", "rho", "sigma"])

# iterate thru our settings
for popsize in popsizes:
    for mutation in mutations:
        for recombination in recombinations:
            
            # iterate through our seeds
            for seed in range(seeds):

                # set our seed
                np.random.seed(seed)

                # define the lorenz system to play nicely with solve_ivp
                def lorenz(t, y, beta, rho, sigma):
                    X, Y, Z = tuple(y)
                    dydt = np.array([sigma*(Y-X), X*(rho-Z)-Y, X*Y - beta*Z])
                    return dydt

                # generate our in-interval data (unnoised) + truncate both as necessary
                ts_interval = np.linspace(0.0, TM, math.ceil((TM - 0.0) * 250 + 1))
                X_ts_interval = solve_ivp(fun=lorenz, t_span=(0.0, TM), y0=y0, 
                                          t_eval=ts_interval, 
                                          args=(beta_true, rho_true, sigma_true)).y.T
                X_ts_interval = X_ts_interval[ts_interval >= TT]
                ts_interval = ts_interval[ts_interval >= TT]

                # compute our noise levels
                noise_levels = alpha*(X_ts_interval.max(axis=0)\
                                      - X_ts_interval.min(axis=0))

                # generate in-sample training data (noised)
                ts_sample = np.linspace(0.0, TM, math.ceil(dobs*(TM - 0.0)))
                X_ts_sample_unnoised = solve_ivp(fun=lorenz, t_span=(0.0, TM), y0=y0, 
                                                 t_eval=ts_sample, 
                                                 args=(beta_true, 
                                                       rho_true, sigma_true)).y.T
                X_ts_sample_noised = X_ts_sample_unnoised +\
                np.random.normal(loc=0.0, scale=noise_levels, 
                                 size=(ts_sample.shape[0], 3))
                X_ts_sample_noised = X_ts_sample_noised[ts_sample >= TT]
                ts_sample = ts_sample[ts_sample >= TT]

                # define our MSE loss function based on the in-sample training data
                # in contrast to PSO, only need to write for ONE INDIVIDUAL!
                def objective_func(thetas):

                    # immediately transpose to facilitate for-looping
                    thetas = thetas.T.reshape(-1, 3)

                    # helper function to compute MSE for one particle
                    def one_particle_obj(theta):

                        # unpack into beta, rho, sigma
                        beta, rho, sigma = tuple(theta)

                        # compute our inferred trajectories + truncate to approp. range
                        Xhat = solve_ivp(fun=lorenz, t_span=(0.0, TM), y0=y0, 
                                         t_eval=ts_sample, args=(beta, rho, sigma)).y.T
                        Xhat = Xhat[ts_sample >= TT]

                        # compute + return our L2 error
                        return ((X_ts_sample_noised - Xhat) ** 2).sum(axis=1).mean()

                    # return a column vector
                    return np.array([one_particle_obj(theta)\
                                     for theta in thetas])


                # instantiate + solve our differential evolution algorithm
                # MUST HAVE REAL-VALUED BOUNDS, CANNOT USE NP.INF!
                res = differential_evolution(func=objective_func, bounds=[(0.0, 30.0)]*3,
                                             maxiter=1000, 
                                             popsize=popsize, mutation=mutation, 
                                             recombination=recombination, seed=seed, vectorized=True, disp=True)

                # extract out our objective function value + theta_hat
                cost, theta_hat = res.fun, res.x

                # create the entry for this run into our bigger dataframe
                row = [popsize, mutation[0], mutation[1], recombination, seed, 
                       theta_hat[0], theta_hat[1], theta_hat[2]]

                # add to our dataframe
                df.loc[len(df.index)] = row

                # save our file PREEMPTIVELY AS A FORM OF CHECKPOINTING!
                df.to_csv(f"results/{regime}/{FNAME}", index=False)