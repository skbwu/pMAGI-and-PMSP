import os
import time
import datetime
import math
import numpy as np
import tensorflow as tf
from decimal import Decimal
from scipy import interpolate
from src import DenseBlock

# main class of their PINN model.
class LorenzPINN(tf.keras.Model):
    
    # constructor -- just inherits keras Model.
    def __init__(self, layers=4, layer_width=20, bn=False, log_opt=False, lr=1e-2, lmbda=10.0):
        
        # initialize
        super(LorenzPINN, self).__init__()

        # sigma, rho, beta in log-scale.
        self.c1 = tf.Variable(1.0, trainable=True, name="log-sigma")
        self.c2 = tf.Variable(1.0, trainable=True, name="log-rho")
        self.c3 = tf.Variable(1.0, trainable=True, name="log-beta")

        # builds their fully-connected layers
        self.NN = DenseBlock(layers, layer_width, bn)
        
        # internal logging of no. of epochs trained
        self.epochs = 0
        
        # are we doing optimization in log-scale or not?
        self.log_opt = log_opt
        
        # what's our optimizer?
        self.optimizer = tf.optimizers.Adam(learning_rate=lr)

        # lmbda is the tradeoff between reconstruction error and physics error
        self.lmbda = lmbda # we'll try 0.1, 1.0, 100.0, 1000.0
        
        # for storing our in-sample value at time TM after the first 60K
        self.TM_vals = np.array([0.0, 0.0, 0.0])
    

    # this is their forward-pass function
    # @tf.function
    def call(self, t, TT, TM, TFC, is_forecasting):
        
        # get our rho, sigma, beta that we optimize in log-space
        theta_hat = self.predict()
        
        # let's get our x, y, z on the in-sample t that's passed in (for use on in-sample recons. error)
        # AS OF 2/20/2024:
        if is_forecasting:
            [x, y, z] = self.NN(t)
            x = (x * ((t - TM) / TFC)) + (self.TM_vals[0] * ((TM + TFC - t) / TFC))
            y = (y * ((t - TM) / TFC)) + (self.TM_vals[1] * ((TM + TFC - t) / TFC))
            z = (z * ((t - TM) / TFC)) + (self.TM_vals[2] * ((TM + TFC - t) / TFC))
        
        # we know we are working with in-sample!
        else:
            [x, y, z] = self.NN(t)
        
        # if we are doing the forecasting phase, our t_physics is ONLY the FORECASTING REGION!
        if is_forecasting:
            
            # AS OF 2/20/2024: we will ALLOW FOR CLOSED INTERVALS ON BOTH SIDES!
            t_physics = tf.convert_to_tensor(np.linspace(start=TM, stop=TM+TFC, 
                                                         num=math.ceil(TFC*40)+1).reshape(-1, 1).astype(np.float32))
        
        # if we're still in the in-sample training phase, our t_physics is IN-SAMPLE ONLY! 40x evenly-spaced pts / unit time.
        else:
            t_physics = tf.convert_to_tensor(np.linspace(start=TT, stop=TM, 
                                                         num=math.ceil((TM-TT)*40)+1).reshape(-1, 1))
        
        # this is the critical component where we must use automatic differentiation!
        with tf.GradientTape(persistent=True) as g:
            
            # we need to watch the t_physics vector as input.
            g.watch(t_physics)
            
            # AS OF 2/20/2024: THIS IS TO HARD-CODE THE BOUNDARY CONSTRAINTS during the forecasting MODE!
            if is_forecasting:
                
                # hard-enforcement of boundary values!
                [x_physics, y_physics, z_physics] = self.NN(t_physics)
                x_physics = (x_physics * ((t_physics - TM) / TFC)) + (self.TM_vals[0] * ((TM + TFC - t_physics) / TFC))
                y_physics = (y_physics * ((t_physics - TM) / TFC)) + (self.TM_vals[1] * ((TM + TFC - t_physics) / TFC))
                z_physics = (z_physics * ((t_physics - TM) / TFC)) + (self.TM_vals[2] * ((TM + TFC - t_physics) / TFC))
            
            # if not forecasting, proceed as usual
            else:
                
                # just directly compute physics error.
                [x_physics, y_physics, z_physics] = self.NN(t_physics)
                

        # gradients w.r.t. t_physics
        dx_dt = g.gradient(x_physics, t_physics)
        dy_dt = g.gradient(y_physics, t_physics)
        dz_dt = g.gradient(z_physics, t_physics)
        
        # note for Prof. Yang -- this seems odd.
        del g

        # calculation of residual between physics-error of autodiff vs. Lorenz system formulas using in-interval
        fx = tf.cast(dx_dt, tf.float32) - theta_hat[0] * (y_physics - x_physics)
        fy = tf.cast(dy_dt, tf.float32) - x_physics * (theta_hat[1] - z_physics) + y_physics
        fz = tf.cast(dz_dt, tf.float32) - x_physics * y_physics + theta_hat[2] * z_physics

        # note that x, y, z are in-sample, while fx, fy, and fz are in-interval.
        return [x, y, z, fx, fy, fz]
    
    
    # function to encode TM_vals
    def store_TM_vals(self, TM_vals):
        self.TM_vals = TM_vals
    
    # just instantiating our learning rate + optimizer
    def set_lr(self, lr):
        self.optimizer = tf.optimizers.Adam(learning_rate=lr)

    # passes into the __mse function
    def get_loss(self, t, TT, TM, TFC, is_forecasting, u_true):
        return self.__mse(u_true, self(t, TT, TM, TFC, is_forecasting))
    
    # this is for checking the parameter recovery only!
    def get_error(self, true):
        pred = tf.squeeze(self.predict())
        true = tf.convert_to_tensor(true)
        return tf.reduce_sum(tf.abs(pred - true))

    # converts our current conceptions of log-sigma, log-rho, log-beta to exponentiated form if necessary
    def predict(self):
        
        # convert to tensor
        var_tensor_mu = tf.convert_to_tensor([self.c1,self.c2,self.c3])
        
        # exponentiate
        exp_var_tensor_mu = tf.exp(var_tensor_mu) if self.log_opt else var_tensor_mu

        return exp_var_tensor_mu

    # generating trajectory reconstruction / forecasting: assuming either in-sample or forecasting!
    def predict_curves(self, t, TM, TFC, is_forecasting):
        
        # cast to column vector to play safe with dimensions
        t = t.reshape(-1, 1)
        
        # use our "NN2"
        if is_forecasting:
            [x, y, z] = self.NN(t)
            x = (x * ((t - TM) / TFC)) + (self.TM_vals[0] * ((TM + TFC - t) / TFC))
            y = (y * ((t - TM) / TFC)) + (self.TM_vals[1] * ((TM + TFC - t) / TFC))
            z = (z * ((t - TM) / TFC)) + (self.TM_vals[2] * ((TM + TFC - t) / TFC))
        
        # we know we are working with in-sample!
        else:
            [x, y, z] = self.NN(t)
        
        # return our output
        return [x, y, z]
    
    # minimize our loss function which is our physics + reconstruction loss weighted sum.
    def optimize(self, t, TT, TM, TFC, is_forecasting, u_true):
        loss = lambda: self.get_loss(t, TT, TM, TFC, is_forecasting, u_true)
        
        # do our gradient step on all the parameters if not forecasting
        if is_forecasting:
            self.optimizer.minimize(loss=loss, var_list=self.trainable_weights[:-3]) # the last three params are theta!
        else:    
            self.optimizer.minimize(loss=loss, var_list=self.trainable_weights)
    
    # the main function called by our pinn runscript.
    def fit(self, observed_data, TT, TM, TFC, is_forecasting, true_pars, epochs, verbose=False):
        
        # for each epoch ... note that we will ALWAYS have epochs = 1 to stay safe!
        for ep in range(self.epochs+1,self.epochs+epochs+1):
            
            # self.optimize calls self.get_loss, which itself is just a wrapper effectively for self.__mse
            self.optimize(observed_data[0], TT, TM, TFC, is_forecasting,
                          [observed_data[1],observed_data[2],observed_data[3]])
        
        # this is an INTERNAL self-counter! we will use the for-loop in our main runner script.
        self.epochs += epochs

    
    # this COMPUTES OUR WEIGHTED RECONSTRUCTION + PHYSICS ERRORS!
    def __mse(self, u_true, y_pred):

        # squared-loss on the reconstruction loss
        loss_x = tf.reduce_mean(tf.square(y_pred[0] - u_true[0]))
        loss_y = tf.reduce_mean(tf.square(y_pred[1] - u_true[1]))
        loss_z = tf.reduce_mean(tf.square(y_pred[2] - u_true[2]))
        
        # squared loss on the physics-error residuals. It is already reduce mean!
        loss_fx = tf.reduce_mean(tf.square(y_pred[3]))
        loss_fy = tf.reduce_mean(tf.square(y_pred[4]))
        loss_fz = tf.reduce_mean(tf.square(y_pred[5]))
    
        # weighted sum of our reconstruction + physics-error residuals.
        return self.lmbda*(loss_x + loss_y + loss_z) + loss_fx + loss_fy + loss_fz




