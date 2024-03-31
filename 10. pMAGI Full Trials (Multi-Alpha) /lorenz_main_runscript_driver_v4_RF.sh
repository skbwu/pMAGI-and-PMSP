# create "outputs" folder if it doesn't exist
if [ ! -d "outputs" ]; then
    mkdir "outputs"
    echo "Created 'outputs' folder."
fi

# create "errors" folder if it doesn't exist
if [ ! -d "errors" ]; then
    mkdir "errors"
    echo "Created 'errors' folder."
fi

# create "results" folder if it doesn't exist
if [ ! -d "results" ]; then
    mkdir "results"
    echo "Created 'results' folder."
fi

# what are our GLOBAL SETTINGS for density of observations, alpha (noise), and specify noises?
dens_arr=("5" "10" "20" "40")
disc_arr=("3" "2" "1" "0")
alpha_arr=("0.15" "0.015" "0.0015" "0.00015" "0.000015")
SN_arr=("FALSE")

# pilot settings: PL = pilot length, PD = pilot discretization
PL_arr=("0.0" "1.0" "2.0")

# get the length of the dens array
array_length=${#dens_arr[@]}

# iterate through the settings
for ((index = 0; index < array_length; index++)); do
    
    # figure out observation density
    dens=${dens_arr[index]}
    
    # figure out what the corresponding discretization should be
    disc=${disc_arr[index]}

    for alpha in "${alpha_arr[@]}"; do
        for SN in "${SN_arr[@]}"; do
            for PL in "${PL_arr[@]}"; do
                
                # determine what our PD_arr should be, based on whether we are specifying noise or not.
                # if PL == "0.0", PD is irrelevant - just set it to overall discretization.
                if [[ "$PL" == "0.0" ]]; then
                    PD_arr=("$disc")
                
                # else if PL != "0.0" but SN == "TRUE", can also just set it to overall discretization.
                elif [[ "$SN" == "TRUE" ]]; then
                    PD_arr=("$disc")
                
                # else if PL != "0.0" but SN == "FALSE", then PD is very, very relevant because we need to infer noise!
                else
                    PD_arr=("0" "1" "2" "3")
                fi
                
                # now iterate through pilot-discretization and launch our SBATCH jobs.
                for PD in "${PD_arr[@]}"; do

                    # 1. stable (canonical), with intervals of lengths 2.0, 4.0, and 10.0
                    sbatch lorenz_main_runscript_v4_RF.sh 6.00 10.00 5.0 5.0 5.0 0.0 2.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 6.00 10.00 5.0 5.0 5.0 0.0 4.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 6.00 10.00 5.0 5.0 5.0 0.0 6.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 6.00 10.00 5.0 5.0 5.0 0.0 8.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 6.00 10.00 5.0 5.0 5.0 0.0 10.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5

                    # 2. stable (transient chaos), with intervals of lengths 2.0, 4.0, and 10.0
                    sbatch lorenz_main_runscript_v4_RF.sh 23.00 10.00 5.0 5.0 5.0 0.0 2.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 23.00 10.00 5.0 5.0 5.0 0.0 4.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 23.00 10.00 5.0 5.0 5.0 0.0 6.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 23.00 10.00 5.0 5.0 5.0 0.0 8.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 23.00 10.00 5.0 5.0 5.0 0.0 10.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5

                    # 3. chaotic (butterfly), with intervals of lengths 2.0, 4.0, and 10.0
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 5.0 5.0 5.0 0.0 2.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 5.0 5.0 5.0 0.0 4.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 5.0 5.0 5.0 0.0 6.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 5.0 5.0 5.0 0.0 8.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 5.0 5.0 5.0 0.0 10.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    
                    # 4. chaotic (no butterfly), with intervals of lengths 2.0, 4.0, and 10.0
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 2.0 2.0 2.0 0.0 2.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 2.0 2.0 2.0 0.0 4.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 2.0 2.0 2.0 0.0 6.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 2.0 2.0 2.0 0.0 8.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    sbatch lorenz_main_runscript_v4_RF.sh 28.00 10.00 2.0 2.0 2.0 0.0 10.0 $dens $disc $alpha $SN 16001 $PL 4001 $PD
                    sleep 0.5
                    
                done
            done
        done        
    done        
done                