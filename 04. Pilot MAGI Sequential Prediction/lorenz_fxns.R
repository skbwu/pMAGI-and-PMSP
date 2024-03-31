# create the function defining Lorenz ODE system
lorenzmodelODE <- function(theta, x, tvec) {
  
  # unpack the variable-timeseries in our case
  X = x[,1]
  Y = x[,2]
  Z = x[,3]
  
  # unpack the parameters
  beta = theta[1]
  rho = theta[2]
  sigma = theta[3]
  
  # create our output of derivatives (should be same shape as 'x')
  XYZdt = array(data=0, dim=c(nrow(x), ncol(x)))
  XYZdt[,1] = sigma * (Y - X)
  XYZdt[,2] = X * (rho - Z) - Y
  XYZdt[,3] = (X * Y) - (beta * Z)
  
  # return the derivatives
  XYZdt
}


# magi needs the derivatives w.r.t. the variables ...
lorenzmodelDx <- function(theta, x, tvec) {
  
  # unpack the variable-timeseries in our case
  X = x[,1]
  Y = x[,2]
  Z = x[,3]
  
  # unpack the parameters
  beta = theta[1]
  rho = theta[2]
  sigma = theta[3]
  
  # note that ArXiV paper's matrix is a bit misleading. Should be rotated 90 degrees
  resultDx = array(0, c(nrow(x), ncol(x), ncol(x)))
  
  # populate each column (first slice)
  resultDx[,1,1] = -sigma
  resultDx[,2,1] = sigma
  resultDx[,3,1] = 0.0
  
  # second slice
  resultDx[,1,2] = rho - Z
  resultDx[,2,2] = -1.0
  resultDx[,3,2] = -X
  
  # third slice
  resultDx[,1,3] = Y
  resultDx[,2,3] = X
  resultDx[,3,3] = -beta
  
  # return our resultDx
  resultDx
}


lorenzmodelDtheta <- function(theta, x, tvec){
  
  # unpack the variable-timeseries in our case
  X = x[,1]
  Y = x[,2]
  Z = x[,3]
  
  # unpack the parameters
  beta = theta[1]
  rho = theta[2]
  sigma = theta[3]
  
  # note that ArXiV paper's matrix is a bit misleading. Should be rotated 90 degrees
  resultDtheta = array(0, c(nrow(x), length(theta), ncol(x)))
  
  # populate each column (first slice)
  resultDtheta[,1,1] = 0.0
  resultDtheta[,2,1] = 0.0
  resultDtheta[,3,1] = Y - X
  
  # second slice
  resultDtheta[,1,2] = 0.0
  resultDtheta[,2,2] = X
  resultDtheta[,3,2] = 0.0
  
  # third slice
  resultDtheta[,1,3] = -Z
  resultDtheta[,2,3] = 0.0
  resultDtheta[,3,3] = 0.0
  
  # return our resultDtheta
  resultDtheta
}


# create our ODE function in a form that's suitable for numerical integration.
modelODE <- function(tvec, state, parameters){
  
  # 'as' forces the input to become a certain class, in this case, a vector.
  # 't' means transpose.
  list(as.vector(lorenzmodelODE(parameters, t(state), tvec)))
}


# function that extracts phi parameters
extract_phi_csv <- function(result, filename) {
  
  # save as .csv
  write.csv(result$phi, filename, row.names = TRUE)
}
               
               
# function that saves a .csv of the inferred parameter values.
extract_estimates_csv <- function(result, filename) {
  
  # get the estimated 0.025 and 0.975 quantiles of the thetas. '2' means cols.
  theta.est <- apply(result$theta, 2,
                     
                     # this is equivalent of lambda function in Python.
                     function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
  
  # names of our parameters
  theta.names <- c("beta", "rho", "sigma")
  
  # set the column names of theta.est to theta.names
  colnames(theta.est) <- theta.names
  
  # set the row names (equivalent of the DataFrame in Pandas)
  rownames(theta.est) <- c("Mean", "2.5%", "97.5%")
  
  # round to three decimal points
  theta.est <- signif(theta.est, 3)
  
  # save as .csv
  write.csv(theta.est, filename, row.names = TRUE)
}

# function that saves a .csv of the inferred phi values.
extract_phis_csv <- function(result, filename) {
  # do the same for the phi's
  phi.est <- result$phi
  phi.names <- c("beta", "rho", "sigma")
  colnames(phi.est) <- phi.names
  write.csv(phi.est, filename, row.names = TRUE)
}


# compile everything into a MAGI payload
lorenzmodel <- list(
  fOde = lorenzmodelODE,
  fOdeDx = lorenzmodelDx,
  fOdeDtheta = lorenzmodelDtheta,
  thetaLowerBound = rep(0.0, 3),
  thetaUpperBound = rep(30.0, 3)
)

