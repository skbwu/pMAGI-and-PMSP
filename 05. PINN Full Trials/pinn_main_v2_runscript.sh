#!/bin/bash
#SBATCH -J PINN2_RF # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p shared # Partition
#SBATCH --mem 16000 # Memory request (64 GB), used to be 64 GB but we can cut down to 16 GB
#SBATCH -t 0-05:00 # (D-HH:MM) # used to be 36 hours, but only doing 1 trial now to make amends.
#SBATCH -o /n/holyscratch01/kou_lab/swu/010223/pinn_sequential_forecasting_inference/outputs/%j.out # Standard output
#SBATCH -e /n/holyscratch01/kou_lab/swu/010223/pinn_sequential_forecasting_inference/errors/%j.err # Standard error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skylerwu@college.harvard.edu
#SBATCH --account=kou_lab

module load cmake/3.25.2-fasrc01
module load gcc/12.2.0-fasrc01
conda run -n afterburner python3 pinn_main_v2.py $1 $2 $3 $4 $5 $6
