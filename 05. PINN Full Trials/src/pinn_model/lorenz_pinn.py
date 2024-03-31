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
        self.c1 = tf.Variable(1.0)
        self.c2 = tf.Variable(1.0)
        self.c3 = tf.Variable(1.0)

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

    # this is their forward-pass function
    # @tf.function
    def call(self, t, TT, TM, deltaT_forecast, fcPhysics):
        
        # get our rho, sigma, beta that we optimize in log-space
        theta_hat = self.predict()
        
        # let's get our x, y, z on the in-sample t that's passed in (for use on in-sample recons. error)
        [x, y, z] = self.NN(t)
        
        # for physics error, we need a finer-grid of observations (40 obs / unit time) -- need to account for pot. forecasting!
        if fcPhysics:
            
            # also append our forecasting observations too!
            t_physics = tf.convert_to_tensor(np.concatenate([np.linspace(start=TT, stop=TM, 
                                                                         num=math.ceil((TM - TT) * 40)),
                                                             np.linspace(start=TM, stop=TM + deltaT_forecast,
                                                                         num=math.ceil(40*deltaT_forecast)+1)[1:]])\
                                             .reshape(-1, 1))
        
        # then, just use our in-sample discretization
        else:
            t_physics = tf.convert_to_tensor(np.linspace(start=TT, stop=TM, 
                                                         num=math.ceil((TM - TT) * 40)).reshape(-1, 1))

        # this is the critical component where we must use automatic differentiation!
        with tf.GradientTape(persistent=True) as g:
            g.watch(t_physics)
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

    
    # just instantiating our learning rate + optimizer
    def set_lr(self, lr):
        self.optimizer = tf.optimizers.Adam(learning_rate=lr)

    # passes into the __mse function
    def get_loss(self, t, TT, TM, deltaT_forecast, fcPhysics, u_true):
        return self.__mse(u_true, self(t, TT, TM, deltaT_forecast, fcPhysics))
    
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

    # generating trajectory reconstruction
    def predict_curves(self, t):
        return self.NN(t)
    
    # minimize our loss function which is our physics + reconstruction loss weighted sum.
    def optimize(self, t, TT, TM, deltaT_forecast, fcPhysics, u_true):
        loss = lambda: self.get_loss(t, TT, TM, deltaT_forecast, fcPhysics, u_true)
        self.optimizer.minimize(loss=loss, var_list=self.trainable_weights)
    
    # the main function called by our pinn runscript.
    def fit(self, observed_data, TT, TM, deltaT_forecast, fcPhysics, true_pars, epochs, verbose=False):
        
        # for each epoch ... note that we will ALWAYS have epochs = 1 to stay safe!
        for ep in range(self.epochs+1,self.epochs+epochs+1):
            
            # self.optimize calls self.get_loss, which itself is just a wrapper effectively for self.__mse
            self.optimize(observed_data[0], TT, TM, deltaT_forecast, fcPhysics,
                          [observed_data[1],observed_data[2],observed_data[3]])
            
            # the below code is simply to tell us how our model is performing over time. Not actually affecting model!
            if ep % 1000 == 0 or ep == 1:
                
                # pred = self.predict()
                loss = self.get_loss(observed_data[0], TT, TM, deltaT_forecast, fcPhysics,
                                     [observed_data[1],observed_data[2],observed_data[3]]) / observed_data[0].shape[0]
                
                error = self.get_error(true_pars)        
                curves = self.predict_curves(observed_data[0])
                
                # purely for console output purposes
                if verbose:
                    print('\n')
                    print(
                        "Epoch: {:5d}, loss: {:.2E}, error: {:3.2f}".format(
                            ep, Decimal(loss.numpy().item()), error.numpy().item()
                        )                    
                    )
                    print(
                        "c1: {:3.2f}, c2: {:3.2f}, c3: {:3.2f}".format(
                            np.exp(self.c1.numpy().item()), np.exp(self.c2.numpy().item()), np.exp(self.c3.numpy().item())
                        )                    
                    )
        
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




