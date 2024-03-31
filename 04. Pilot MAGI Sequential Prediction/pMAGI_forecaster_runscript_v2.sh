#!/bin/bash
#SBATCH -J PMFC2 # A single job name for the array
#SBATCH -n 1 # One core.
#SBATCH -N 1 # One node.
#SBATCH -p sapphire,shared # Partition
#SBATCH --mem 64000 # Memory request (64 GB)
#SBATCH -t 1-06:00 # (D-HH:MM)
#SBATCH -o /n/home11/skbwu/06p2_pMAGI_sequential-forecasting_revised/outputs/PMFC2_%j.out
#SBATCH -e /n/home11/skbwu/06p2_pMAGI_sequential-forecasting_revised/errors/PMFC2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skylerwu@college.harvard.edu

module load cmake/3.25.2-fasrc01
module load gcc/12.2.0-fasrc01
module load R/4.2.2-fasrc01 # load in the R module
export R_LIBS_USER=/n/home11/skbwu/apps/R_422 # tell R where to look for locally installed packages

# call in our scripts - there should only be three settings
Rscript pMAGI_forecaster_v2.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
