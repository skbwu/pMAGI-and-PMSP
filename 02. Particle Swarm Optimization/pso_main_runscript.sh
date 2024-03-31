#!/bin/bash
#SBATCH -J PSO # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p shared # Partition
#SBATCH --mem 16000 # Memory request (16 GB)
#SBATCH -t 0-24:00 # (D-HH:MM)
#SBATCH -o /n/holyscratch01/kou_lab/swu/122323/PSO/outputs/%j.out # Standard output
#SBATCH -e /n/holyscratch01/kou_lab/swu/122323/PSO/errors/%j.err # Standard error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skylerwu@college.harvard.edu

# just directly run the script!
conda run -n afterburner python3 pso_main.py $1 $2 $3