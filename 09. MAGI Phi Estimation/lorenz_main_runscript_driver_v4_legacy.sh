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

# what are our GLOBAL SETTINGS for alpha (noise) and seed?
alpha_arr=("0.15" "0.015" "0.0015" "0.00015" "0.000015")
seed_arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")

# just need to iterate over regime (manually), noise level, and seed
for alpha in "${alpha_arr[@]}"; do
    for seed in "${seed_arr[@]}"; do
    
        # one script call per regime.
        sbatch lorenz_main_runscript_v4_legacy.sh 6.00 5.0 5.0 5.0 $alpha $seed
        sleep 0.5
        sbatch lorenz_main_runscript_v4_legacy.sh 23.00 5.0 5.0 5.0 $alpha $seed
        sleep 0.5
        sbatch lorenz_main_runscript_v4_legacy.sh 28.00 5.0 5.0 5.0 $alpha $seed
        sleep 0.5
        sbatch lorenz_main_runscript_v4_legacy.sh 28.00 2.0 2.0 2.0 $alpha $seed
        sleep 0.5

    done
done