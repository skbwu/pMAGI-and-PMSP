#!/bin/bash
#SBATCH -J LM4_RF # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p sapphire,shared # Partition
#SBATCH --mem 16000 # Memory request (16 GB)
#SBATCH -t 0-12:00 # (D-HH:MM)
#SBATCH -o /n/holyscratch01/kou_lab/swu/camera_ready/01.2_pMAGI+MAGI_full-trials_REFIRE/outputs/LMRF_%j.out # Standard output
#SBATCH -e /n/holyscratch01/kou_lab/swu/camera_ready/01.2_pMAGI+MAGI_full-trials_REFIRE/errors/LMRF_%j.err # Standard error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skylerwu@college.harvard.edu

module load cmake/3.25.2-fasrc01
module load gcc/12.2.0-fasrc01
module load R/4.2.2-fasrc01 # load in the R module
export R_LIBS_USER=/n/home11/skbwu/apps/R_422 # tell R where to look for locally installed packages

# run 10x trials for this given setting
for seed in {0..9}; do

    # call in our scripts - there should only be three settings
    Rscript lorenz_main_v4_RF.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} "$seed" # directly run the script!
    
    # make sure the script finishes running before we launch the next one
    wait

done