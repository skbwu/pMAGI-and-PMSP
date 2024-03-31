#!/bin/bash
#SBATCH -J PFC_v2 # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p sapphire,shared # Partition
#SBATCH --mem 16000 # Memory request (16 GB)
#SBATCH -t 0-06:00 # (D-HH:MM)
#SBATCH -o /n/holyscratch01/kou_lab/swu/camera_ready/07p3_PINN_sequential-forecasting_revised/outputs/%j.out # Standard output
#SBATCH -e /n/holyscratch01/kou_lab/swu/camera_ready/07p3_PINN_sequential-forecasting_revised/errors/%j.err # Standard error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skylerwu@college.harvard.edu
#SBATCH --account=kou_lab

module load cmake/3.25.2-fasrc01
module load gcc/12.2.0-fasrc01
conda run -n afterburner python3 pinn_forecaster_main.py $1 $2 $3 $4 $5 $6
