# Contents and Overview
This repository contains the source code for reproducing the main experiments for the undergraduate senior thesis 
*Extracting Signal out of Chaos: Advancements on MAGI for Bayesian Analysis of Dynamical Systems* written by Skyler Wu and advised by Professor Samuel Kou of the Department of Statistics at Harvard University.

Each folder in this repository corresponds to one major experiment / model class in *Extracting Signal out of Chaos*. 
The contents of this repository reproduce the raw model outputs (e.g., HMC samples, mean trajectories, point estimates of parameters, etc.) 
of all models discussed in *Extracting Signal out of Chaos*, including MAGI, Pilot MAGI (pMAGI), Pilot MAGI Sequential Prediction (PMSP), 
Particle Swarm Evolution (PSO), Differential Evolution (DE), and Physics-Informed Neural Networks PINNs). This repository also includes the datasets used to benchmark competitor methods against pMAGI and PMSP, 
including both the unnoised ground truth and the noised observations.

# Running Experiments
Many of the experiments performed for this senior thesis were run on the FASRC cluster supported by the FAS Division of Science Research Computing Group at Harvard University. 
As such, most experiment folders in this repository were designed to run on a SLURM-based cluster. Almost every folder contains at least one file of the form `(...)_runscript_driver.sh`. 
Within an appropriate SLURM cluster environment (adjusting partitions, fairshare accounts, etc. accordingly, of course), one can call `bash (...)_runscript_driver.sh` to send a set of jobs
 containing every hyperparameter/data variant of an experiment/model to the cluster.

 # Questions?
 Any questions regarding this repository can be directed to [skylerwu@college.harvard.edu](skylerwu@college.harvard.edu).
