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

# iterate through the settings
for ((a = 0; a < 4; a++)); do
    for ((b = 0; b < 4; b++)); do
        for ((c = 0; c < 2; c++)); do
            for ((d = 0; d < 5; d++)); do
                for ((e = 0; e < 5; e++)); do
                    for ((f = 0; f < 5; f++)); do

                        # launch our job
                        sbatch pinn_forecaster_main_conservative_runscript.sh $a $b $c $d $e $f
                        sleep 0.5
                    
                    done
                done
            done
        done        
    done        
done                