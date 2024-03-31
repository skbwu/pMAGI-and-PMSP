import numpy as np
import tensorflow as tf
import pandas as pd
from src import LorenzPINN # that's where their model is
from scipy.integrate import solve_ivp
import sys, copy, os, shutil, time, math


'''
let's try PINN on the following settings
1. 10x randomized trials per setting.
2. 2 stable + 2 chaotic regimes.
3. dens_obs: ONLY 5, 10, 20, 40
4. noise levels: ONLY 1.5e-1
5. interval lengths of: 2, 4, 6, 8, 10, starting from TT=0.0 only.

*evaluate on a grid of tmax*250+1 points.
'''

# for purposes of notebook only!
# sys.argv = [sys.argv[0], "2", "3", "0", "4", "2", "0"]

# what are our possible regime settings?
regime_descs = ["stable_canonical", "stable_transient-chaos", "chaotic_butterfly", "chaotic_no-butterfly"]
theta_vals = [(8/3, 6.0, 10.0), (8/3, 23.0, 10.0), (8/3, 28.0, 10.0), (8/3, 28.0, 10.0)]
init_conds = [(5.0, 5.0, 5.0), (5.0, 5.0, 5.0), (5.0, 5.0, 5.0), (2.0, 2.0, 2.0)]
lmbdas = [0.1, 1.0, 10.0, 100.0, 1000.0]
fcPhysics = [True, False] # are we incorporating forecasting range into physics-error computation?

# data availability + scope settings.
d_obs_list = [5, 10, 20, 40]
alphas = [1.5e-1]
intervals = [(0.0, 2.0), (0.0, 4.0), (0.0, 6.0), (0.0, 8.0), (0.0, 10.0)]

# regime-specific command line arguments
regime = regime_descs[int(sys.argv[1])]
beta, rho, sigma = theta_vals[int(sys.argv[1])]
y0 = np.array(init_conds[int(sys.argv[1])])

# data noisyness + availability settings / interval-size command line args
d_obs = d_obs_list[int(sys.argv[2])]
alpha = alphas[int(sys.argv[3])]
TT, TM = intervals[int(sys.argv[4])]

# PINN: also need to specify our lmbda (reconstruction vs. physics error tradeoff), and fcPhysics
lmbda = lmbdas[int(sys.argv[5])]
fcPhysics = False # just hardcode it! fcPhysics[int(sys.argv[6])]

# MOST IMPORTANT SETTING FOR THIS FORECASTING EXPERIMENT!
deltaT_forecast = 5.0


# fixed number of epochs + trials
n_epochs = 60000
seeds = 10

# iterate through our seeds (MEMORY ISSUES ONLY AROSE ON THE LAST SEED!)
for seed in [9]: # range(seeds):
    
    # name our directory for storing our results + create it
    dir_name = f"results/{regime}/dobs={d_obs}_alpha={alpha}_lmbda={lmbda}_FCP={fcPhysics}_TT={TT}_TM={TM}_TFC={deltaT_forecast}_epochs={n_epochs}/seed={seed}"
    os.makedirs(dir_name, exist_ok=True)
    
    #### DATA-GENERATION WITH NOISE! SET A SEED!
    np.random.seed(seed)

    # Define the Lorenz system as a function of (x, y, z), t and rho
    def lorenz(t, y, beta, rho, sigma):
        X, Y, Z = tuple(y)
        dydt = np.array([sigma*(Y-X), X*(rho-Z)-Y, X*Y - beta*Z])
        return dydt

    # generate our in-interval data (unnoised) + truncate both as necessary
    ts_interval = np.linspace(0.0, TM, math.ceil((TM - 0.0) * 250 + 1))
    X_ts_interval = solve_ivp(fun=lorenz, t_span=(0.0, TM), y0=y0, 
                              t_eval=ts_interval, args=(beta, rho, sigma)).y.T
    X_ts_interval = X_ts_interval[ts_interval >= TT]
    ts_interval = ts_interval[ts_interval >= TT]

    # compute our noise levels + generate our data
    noise_levels = alpha*(X_ts_interval.max(axis=0) - X_ts_interval.min(axis=0))

    # generate our in-sample training data (noised) + make sure to add the noise
    ts_sample = np.linspace(0.0, TM, math.ceil(d_obs*(TM - 0.0)))
    X_ts_sample_unnoised = solve_ivp(fun=lorenz, t_span=(0.0, TM), y0=y0, 
                                     t_eval=ts_sample, args=(beta, rho, sigma)).y.T    
    X_ts_sample_noised = X_ts_sample_unnoised + np.random.normal(loc=0.0, scale=noise_levels, 
                                                                 size=(ts_sample.shape[0], 3))
    X_ts_sample_noised = X_ts_sample_noised[ts_sample >= TT]
    ts_sample = ts_sample[ts_sample >= TT]

    # create our lorenz_data structure that the PINN needs (NOISED)
    lorenz_data = [ts_sample.reshape(-1, 1), 
                   X_ts_sample_noised[:,0:1], 
                   X_ts_sample_noised[:,1:2], 
                   X_ts_sample_noised[:,2:3]]

    # define what our true parameters sigma, rho, and beta are 
    pars = [sigma, rho, beta]
    
    # solve the entire system including the to-be-forecasted interval: SAME TEST SET FOR ALL SETTINGS!
    # 40 dobs / unit time.
    ts_forecast = np.linspace(TM, TM + deltaT_forecast, math.ceil(40*deltaT_forecast)+1)[1:]
    ts_forecast_combined = np.concatenate([ts_sample, ts_forecast])
    X_ts_forecast_combined_unnoised = solve_ivp(fun=lorenz, t_span=(0.0, TM + deltaT_forecast), y0=y0, 
                                                t_eval=ts_forecast_combined, args=(beta, rho, sigma)).y.T  
    X_ts_forecast_combined_unnoised = X_ts_forecast_combined_unnoised[ts_forecast_combined >= TT]
    ts_forecast_combined = ts_forecast_combined[ts_forecast_combined >= TT]
    X_ts_forecast = X_ts_forecast_combined_unnoised[ts_forecast_combined > TM] # just the fcst region!
    
    '''
    For reproducibility, let's save our data:
    a) in-interval unnoised data
    b) in-sample noised data
    c) in-sample unnoised data
    d) settings of parameters, noise levels, density of observations, etc.
    '''
    
    # a) in-interval unnoised data
    pd.DataFrame(np.hstack([ts_interval.reshape(-1, 1), X_ts_interval]), 
                 columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/y_interval.csv", index=False)
    
    # b) in-sample unnoised data
    pd.DataFrame(np.hstack([ts_sample.reshape(-1, 1), X_ts_sample_noised]), 
                 columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/y_samp_noised.csv", index=False)
    
    # c) in-sample unnoised data
    pd.DataFrame(np.hstack([ts_sample.reshape(-1, 1), X_ts_sample_unnoised]), 
                 columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/y_samp_unnoised.csv", index=False)
    
    
    # d) in-sample + forecast-range's unnoised data
    pd.DataFrame(np.hstack([ts_forecast_combined.reshape(-1, 1), X_ts_forecast_combined_unnoised]), 
                 columns=["t", "X", "Y", "Z"])\
    .to_csv(dir_name + "/y_forecast_combined_unnoised.csv", index=False)
    
    # e) save our settings
    pd.DataFrame(data=np.array([regime, beta, rho, sigma] + list(y0)\
                               + [d_obs, alpha, TT, TM, deltaT_forecast, fcPhysics]).reshape(1, -1),
                 columns=["regime", "beta", "rho", "sigma", "X0", \
                          "Y0", "Z0", "d_obs", "alpha", "TT", "TM", "TFC", "FCP"])\
    .to_csv(dir_name + "/settings.csv", index=False)
    
    ##### PINN SIMULATIONS
    
    # create a dataframe to store our results
    logs = pd.DataFrame(data=None, columns=["epoch", "time", "beta", "rho", "sigma", 
                                            "L1_X_in-sample_noised", "L1_Y_in-sample_noised", "L1_Z_in-sample_noised", 
                                            "L2_X_in-sample_noised", "L2_Y_in-sample_noised", "L2_Z_in-sample_noised",
                                            "L1_X_in-sample_unnoised", "L1_Y_in-sample_unnoised", "L1_Z_in-sample_unnoised", 
                                            "L2_X_in-sample_unnoised", "L2_Y_in-sample_unnoised", "L2_Z_in-sample_unnoised",
                                            "L1_X_in-interval_unnoised", "L1_Y_in-interval_unnoised", "L1_Z_in-interval_unnoised", 
                                            "L2_X_in-interval_unnoised", "L2_Y_in-interval_unnoised", "L2_Z_in-interval_unnoised",
                                            "L1_X_forecast_unnoised", "L1_Y_forecast_unnoised", "L1_Z_forecast_unnoised", 
                                            "L2_X_forecast_unnoised", "L2_Y_forecast_unnoised", "L2_Z_forecast_unnoised",])

    # set a seed for reproducibility
    tf.random.set_seed(seed)

    # they are doing batch-norm and log-transforming all of their beta, rho, and sigma (showcased settings)
    pinn = LorenzPINN(bn=True, log_opt=True, lr=1e-2, layers=3, 
                      layer_width=32, lmbda=lmbda)
    
    # train the PINN for 60K epochs.
    for i in range(n_epochs):
        
        # start our timer
        start = time.time()
        
        # one optimization step
        pinn.fit(lorenz_data, TT, TM, deltaT_forecast, fcPhysics, 
                 pars, 1, verbose=False)

        # end our timer + record time_elapsed
        end = time.time()
        time_elapsed = end - start

        '''
        Metrics of interest: 
        a) epoch, 
        b) time per epoch, 
        c) thetas
        d) L1 in-sample reconstruction error (post-hoc scaled) - noised
        e) L2 in-sample reconstruction error (post-hoc scaled) - noised
        f) L1 in-sample reconstruction error (post-hoc scaled) - unnoised
        g) L2 in-sample reconstruction error (post-hoc scaled) - unnoised
        h) L1 in-interval reconstruction error (post-hoc scaled) - unnoised
        i) L2 in-interval reconstruction error (post-hoc scaled) - unnoised
        '''
        # just record every 100 and the last one
        if (i % 100 == 0) or (i == n_epochs - 1):
        
            # start our row with a) and b)
            row = [i+1, time_elapsed]

            # c) extract our thetas + add to our row
            beta_hat, rho_hat, sigma_hat = np.exp(pinn.c3.numpy()), np.exp(pinn.c2.numpy()), np.exp(pinn.c1.numpy())
            row += [beta_hat, rho_hat, sigma_hat]

            # d) + e) compute in-sample trajectories + compute errors w.r.t NOISED IN-SAMPLE
            curves_sample = pinn.predict_curves(lorenz_data[0])
            Xhat_ts_sample = np.hstack([curve.numpy() for curve in curves_sample])
            row += list(np.abs(X_ts_sample_noised - Xhat_ts_sample).mean(axis=0) / np.abs(X_ts_sample_noised).mean(axis=0))
            row += list(np.sqrt(((X_ts_sample_noised - Xhat_ts_sample) ** 2).mean(axis=0)) / np.sqrt((X_ts_sample_noised ** 2).mean(axis=0)))

            # f) + g) using same COMPUTED IN-SAMPLE TRAJECTORIES, compute errors w.r.t UNNOISED IN-SAMPLE
            row += list(np.abs(X_ts_sample_unnoised - Xhat_ts_sample).mean(axis=0) / np.abs(X_ts_sample_unnoised).mean(axis=0))
            row += list(np.sqrt(((X_ts_sample_unnoised - Xhat_ts_sample) ** 2).mean(axis=0)) / np.sqrt((X_ts_sample_unnoised ** 2).mean(axis=0)))

            # h) + i) compute IN-INTERVAL TRAJECTORIES + compute errors w.r.t UNNOISED IN-INTERVAL
            curves_interval = pinn.predict_curves(ts_interval)
            Xhat_ts_interval = np.hstack([curve.numpy() for curve in curves_interval])
            row += list(np.abs(X_ts_interval - Xhat_ts_interval).mean(axis=0) / np.abs(X_ts_interval).mean(axis=0))
            row += list(np.sqrt(((X_ts_interval - Xhat_ts_interval) ** 2).mean(axis=0)) / np.sqrt((X_ts_interval ** 2).mean(axis=0)))
            
            # j) + k) compute FORECASTED TRAJECTORIES + compute errors w.r.t FORECASTED range
            curves_forecast = pinn.predict_curves(ts_forecast)
            Xhat_ts_forecast = np.hstack([curve.numpy() for curve in curves_forecast])
            row += list(np.abs(X_ts_forecast - Xhat_ts_forecast).mean(axis=0) / np.abs(X_ts_forecast).mean(axis=0))
            row += list(np.sqrt(((X_ts_forecast - Xhat_ts_forecast) ** 2).mean(axis=0)) / np.sqrt((X_ts_forecast ** 2).mean(axis=0)))
            
            # add our row to our logs
            logs.loc[len(logs.index)] = row  
    
    # after all our epochs, save our logs
    logs.to_csv(dir_name + "/logs.csv", index=False)

    # save our final in-sample predictions
    pd.DataFrame(np.hstack([ts_sample.reshape(-1, 1), Xhat_ts_sample]), 
                 columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/yhat_sample.csv", index=False)

    # save our final in-interval predictions
    pd.DataFrame(np.hstack([ts_interval.reshape(-1, 1), Xhat_ts_interval]), 
                 columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/yhat_interval.csv", index=False)
    
    # save our final FORECASTED predictions
    pd.DataFrame(np.hstack([ts_forecast.reshape(-1, 1), Xhat_ts_forecast]), 
                 columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/yhat_forecast.csv", index=False)