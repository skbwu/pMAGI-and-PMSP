#!/bin/bash
#SBATCH -J PHIB # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Number of nodes
#SBATCH -p shared,sapphire # Partition
#SBATCH --mem 8000 # Memory request (16 GB)
#SBATCH -t 0-06:00 # (D-HH:MM)
#SBATCH -o /n/holyscratch01/kou_lab/swu/camera_ready/08_MAGI_phi-estimation/outputs/PHIB_%j.out # Standard output
#SBATCH -e /n/holyscratch01/kou_lab/swu/camera_ready/08_MAGI_phi-estimation/errors/PHIB_%j.err # Standard error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skylerwu@college.harvard.edu

module load cmake/3.25.2-fasrc01
module load gcc/12.2.0-fasrc01
module load R/4.2.2-fasrc01 # load in the R module
export R_LIBS_USER=/n/home11/skbwu/apps/R_422 # tell R where to look for locally installed packages


# array of interval lengths, density of observations, and discretizations.
TM_arr=("0.5" "1.0" "1.5" "2.0" "2.5" "3.0" "3.5" "4.0" "4.5" "5.0" "5.5" "6.0" "6.5" "7.0" "7.5" "8.0" "8.5" "9.0" "9.5" "10.0")
dens_arr=("5" "10" "20" "40")
disc_arr=("3" "2" "1" "0")

# get the length of the dens array
array_length=${#dens_arr[@]}

# need to iterate over interval lengths, and then observation density & corresponding discretization

# iterate through the settings
for ((index = 0; index < array_length; index++)); do
    
    # figure out observation density
    dens=${dens_arr[index]}
    
    # figure out what the corresponding discretization should be
    disc=${disc_arr[index]}

    # now iterate through the TM interval sizes.
    for TM in "${TM_arr[@]}"; do
    
        # call in our scripts - there should only be three settings
        Rscript lorenz_main_v4_legacy.R $1 $2 $3 $4 $TM $dens $disc $5 $6
        
        # make sure the script finishes running before we launch the next one
        wait

    done 
done     
