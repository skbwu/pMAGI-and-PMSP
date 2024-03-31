import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math, sys, copy, os, shutil
import pyswarms as ps
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

# PSO-specific settings: there was one 2007 paper that did c1, c2, w = 2, 2, (0.9 to 0.4)
# 120 particles when doing all three params, 40 particles when doing 2 params, 20 particles with 1 param
# https://www.sciencedirect.com/science/article/abs/pii/S0960077906003067
# also restricted sigma \in (9, 11), beta \in (2, 3), rho \in (20, 30)
dimensions, c1, c2, n_particles = 3, 2.0, 2.0, 120
ws = np.arange(0.9, 0.31, -0.1) # 6 values

# create a directory to store our results + subdirectory for this regime.
if "results" not in os.listdir(): os.mkdir("results")
if regime not in os.listdir("results"): os.mkdir(f"results/{regime}")
    
# extract our true beta, rho, sigma
beta_true, rho_true, sigma_true = theta_true

# create the name of the .csv file that we are going to use for this set of trials.
FNAME = f"TM={TM}_alpha={alpha}_dobs={dobs}.csv"

# create our dataframe to store our results for this Lorenz sim, cols = PSO settings.
df = pd.DataFrame(data=None, columns=["n_particles", "c1", "c2", "w", 
                                      "seed", "beta", "rho", "sigma"])

# iterate through our specific PSO settings: just w!
for w in ws:

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
        # NOTE THAT THIS IS COMPUTED FOR EACH PARTICLE IN THE SWARM!
        def objective_func(theta_particles):

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
                             for theta in theta_particles])


        # instantiate our PSO
        options = {'c1': c1, 'c2': c2, 'w': w}
        optimizer = ps.single.GlobalBestPSO(n_particles=n_particles, 
                                            dimensions=dimensions, 
                                            options=options,
                                            bounds=(np.zeros(3,), np.ones(3,)*30.0))

        # run our PSO + extract results: comment verbose=True to see outputs / status!
        # IT WILL CRASH IF DONT SPECIFY BOUNDS AS SHOWN ABOVE!
        cost, theta_hat = optimizer.optimize(objective_func=objective_func, 
                                             iters=1000, 
                                             verbose=False)

        # create the entry for this run into our bigger dataframe
        row = [n_particles, c1, c2, w, seed, 
               theta_hat[0], theta_hat[1], theta_hat[2]]

        # add to our dataframe
        df.loc[len(df.index)] = row
        
        # save our file IMMEDIATELY for checkpointing purposes!
        df.to_csv(f"results/{regime}/{FNAME}", index=False)