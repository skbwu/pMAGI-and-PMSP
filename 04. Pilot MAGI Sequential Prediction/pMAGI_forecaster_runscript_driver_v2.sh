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

# global settings to iterate over?
TM_arr=("2.0" "10.0" "6.0" "8.0" "4.0")
SN_arr=("FALSE" "TRUE")
DOBS_arr=("40" "20" "10")
ALPHA_arr=("0.000015" "0.15")
STEPSIZE_arr=("0.5" "1.0" "5.0")
SEED_arr=("0")
PM_arr=("1" "5" "4" "3" "2")
DISC_arr=("2" "1")


# blast thru all of our settings!
for SEED in "${SEED_arr[@]}"; do
    for TM in "${TM_arr[@]}"; do
        for SN in "${SN_arr[@]}"; do
            for DOBS in "${DOBS_arr[@]}"; do
                for ALPHA in "${ALPHA_arr[@]}"; do
                    for STEPSIZE in "${STEPSIZE_arr[@]}"; do
                        for PM in "${PM_arr[@]}"; do
                            for DISC in "${DISC_arr[@]}"; do

                                # launch each regime manually
                                sbatch pMAGI_forecaster_runscript_v2.sh 6.00 5.0 $TM $SN $DOBS $DISC $ALPHA $STEPSIZE $PM $SEED
                                sleep 0.5
                                sbatch pMAGI_forecaster_runscript_v2.sh 23.00 5.0 $TM $SN $DOBS $DISC $ALPHA $STEPSIZE $PM $SEED
                                sleep 0.5
                                sbatch pMAGI_forecaster_runscript_v2.sh 28.00 5.0 $TM $SN $DOBS $DISC $ALPHA $STEPSIZE $PM $SEED
                                sleep 0.5
                                sbatch pMAGI_forecaster_runscript_v2.sh 28.00 2.0 $TM $SN $DOBS $DISC $ALPHA $STEPSIZE $PM $SEED
                                sleep 0.5
                            
                            done
                        done
                    done
                done        
            done  
        done
    done     
done
