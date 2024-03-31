import numpy as np
import tensorflow as tf
import pandas as pd
from src import LorenzPINN
import sys, copy, os, shutil, time, math

# what are our possible regime settings?
regime_descs = ["stable_canonical", "stable_transient-chaos", "chaotic_butterfly", "chaotic_no-butterfly"]
theta_vals = [(8/3, 6.0, 10.0), (8/3, 23.0, 10.0), (8/3, 28.0, 10.0), (8/3, 28.0, 10.0)]
init_conds = [(5.0, 5.0, 5.0), (5.0, 5.0, 5.0), (5.0, 5.0, 5.0), (2.0, 2.0, 2.0)]
lmbdas = [0.1, 1.0, 10.0, 100.0, 1000.0]

# data availability + scope settings.
d_obs_list = [5, 10, 20, 40]
alphas = [0.15, 1.5e-5]
intervals = [(0.0, 2.0), (0.0, 4.0), (0.0, 6.0), (0.0, 8.0), (0.0, 10.0)]

# regime-specific command line arguments
regime = regime_descs[int(sys.argv[1])]
beta, rho, sigma = theta_vals[int(sys.argv[1])]
y0 = np.array(init_conds[int(sys.argv[1])])

# data noisyness + availability settings / interval-size command line args
d_obs = d_obs_list[int(sys.argv[2])]
alpha = alphas[int(sys.argv[3])]
TT, TM = intervals[int(sys.argv[4])]

# PINN: also need to specify our lmbda (reconstruction vs. physics error tradeoff)
lmbda = lmbdas[int(sys.argv[5])]

# MOST IMPORTANT SETTING FOR THIS FORECASTING EXPERIMENT! how many steps are we forecasting into the future?
TFC = 5.0

# fixed number of epochs x2 for in-sample + post-freeze (let's only do one seed because limited compute)
n_epochs = 60000
seed = int(sys.argv[6])

# name our directory for storing our results
dir_name = f"results/{regime}/dobs={d_obs}_alpha={alpha}_lmbda={lmbda}" +\
f"_TT={TT}_TM={TM}_TFC={TFC}_epochs={n_epochs}/seed={seed}"

# create our directory
os.makedirs(dir_name, exist_ok=True)


# load in our unnoised + noised data for this run
data_dir = f"rho={int(rho)}_sigma={int(sigma)}_X0={int(y0[0])}_Y0={int(y0[1])}_Z0={int(y0[2])}" +\
f"_TT=0_TM=10_TFC=5_DOBS=40_alpha={alpha}_seed={seed}"
y_samp_forecast_noised = pd.read_csv(f"forecast_data/{data_dir}/y_samp+forecast_noised.csv")

# implicitly enforce our d_obs
X_ts_sample_noised = y_samp_forecast_noised.query(f"t <= {TM}")
X_ts_sample_noised = X_ts_sample_noised[X_ts_sample_noised.index % (40/d_obs) == 0.0]

# extract our time vector and our component observations
ts_sample = X_ts_sample_noised.t.values
X_ts_sample_noised = X_ts_sample_noised[["X", "Y", "Z"]].to_numpy()

# create our lorenz_data structure that the PINN needs (USING NOISED OBSERVATIONS)
lorenz_data = [ts_sample.reshape(-1, 1), 
               X_ts_sample_noised[:,0:1], 
               X_ts_sample_noised[:,1:2], 
               X_ts_sample_noised[:,2:3]]

# define what our true parameters sigma, rho, and beta are 
pars = [sigma, rho, beta]

# also determine what our forecasting ts are. FIXED 2/15/2024
ts_forecast = y_samp_forecast_noised.query(f"t > {TM}").query(f"t <= {TM + TFC}").t.values

##### CONSTRUCT OUR LOGGING + BUILD OUR PINN

# set a seed for reproducibility
tf.random.set_seed(seed)

# they are doing batch-norm and log-transforming all of their beta, rho, and sigma (showcased settings)
pinn = LorenzPINN(bn=True, log_opt=True, lr=1e-2, layers=3, 
                  layer_width=32, lmbda=lmbda)


##### TRAINING THE PINN's theta + weight parameters on in-sample

# create a dataframe to store our results
logs = pd.DataFrame(data=None, columns=["epoch", "time", "beta", "rho", "sigma"])

# train the PINN for another 60K epochs ON IN-SAMPLE OBSERVATIONS ONLY!
for i in range(n_epochs):

    # start our timer
    start = time.time()

    # one optimization step (we're not forecasting!)
    pinn.fit(observed_data=lorenz_data, 
             TT=TT, TM=TM, TFC=TFC, is_forecasting=False, 
             true_pars=pars, epochs=1, verbose=False)

    # end our timer + record time_elapsed
    end = time.time()
    time_elapsed = end - start

    '''
    Metrics of interest (keep it simple!)
    a) epoch, 
    b) time per epoch, 
    c) thetas
    '''
    # just record every 100 and the last one
    if (i % 100 == 0) or (i == n_epochs - 1):

        # start our row with a) and b)
        row = [i+1, time_elapsed]

        # c) extract our thetas + add to our row
        beta_hat, rho_hat, sigma_hat = np.exp(pinn.c3.numpy()), np.exp(pinn.c2.numpy()), np.exp(pinn.c1.numpy())
        row += [beta_hat, rho_hat, sigma_hat]

        # add our row to our logs
        logs.loc[len(logs.index)] = row  

# after all our epochs, save our logs
logs.to_csv(dir_name + "/in-sample_logs.csv", index=False)

# make + save our final in-sample predictions
curves_sample = pinn.predict_curves(lorenz_data[0])
Xhat_ts_sample = np.hstack([curve.numpy() for curve in curves_sample])
pd.DataFrame(np.hstack([ts_sample.reshape(-1, 1), Xhat_ts_sample]), 
             columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/yhat_sample_pre-forecast.csv", index=False)

# save our trained-on-in-sample weights 
pinn.NN.save_weights(dir_name + "/pinn_in-sample/pinn.ckpt")


##### warm-start PINN's weight parameters for forecasting, with theta frozen, 100% physics error.

# DO NOT toggle the lmbda to have zero L2 reconstruction error weight
# pinn.lmbda = 0.0

# create a dataframe to store our results for the forecasting component
logs = pd.DataFrame(data=None, columns=["epoch", "time", "beta", "rho", "sigma"])

# train the PINN for another 60K epochs ONLY USING PHYSICS-ERROR ON THE FORECASTED REGION!
for i in range(n_epochs):

    # start our timer
    start = time.time()

    # one optimization step (YES FORECASTING!)
    # freeze the theta parameters at their in-sample fitted values. already controled by is_forecasting!
    pinn.fit(observed_data=lorenz_data, 
             TT=TT, TM=TM, TFC=TFC, is_forecasting=True, 
             true_pars=pars, epochs=1, verbose=False)

    # end our timer + record time_elapsed
    end = time.time()
    time_elapsed = end - start

    '''
    Metrics of interest (keep it simple!)
    a) epoch, 
    b) time per epoch, 
    c) thetas
    '''
    # just record every 100 and the last one
    if (i % 100 == 0) or (i == n_epochs - 1):

        # start our row with a) and b)
        row = [i+1, time_elapsed]

        # c) extract our thetas + add to our row
        beta_hat = np.exp(pinn.c3.numpy())
        rho_hat = np.exp(pinn.c2.numpy())
        sigma_hat = np.exp(pinn.c1.numpy())
        row += [beta_hat, rho_hat, sigma_hat]

        # add our row to our logs
        logs.loc[len(logs.index)] = row  

# after all our epochs, save our logs
logs.to_csv(dir_name + "/forecast_logs.csv", index=False)

# make + save our final in-sample predictions
curves_sample = pinn.predict_curves(lorenz_data[0])
Xhat_ts_sample = np.hstack([curve.numpy() for curve in curves_sample])
pd.DataFrame(np.hstack([ts_sample.reshape(-1, 1), Xhat_ts_sample]), 
             columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/yhat_sample_post-forecast.csv", index=False)

# save our final FORECASTED predictions
curves_forecast = pinn.predict_curves(ts_forecast)
Xhat_ts_forecast = np.hstack([curve.numpy() for curve in curves_forecast])
pd.DataFrame(np.hstack([ts_forecast.reshape(-1, 1), Xhat_ts_forecast]), 
             columns=["t", "X", "Y", "Z"]).to_csv(dir_name + "/yhat_forecast.csv", index=False)

# save our trained-on-in-sample weights 
pinn.NN.save_weights(dir_name + "/pinn_forecast/pinn.ckpt")