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

# iterate through the settings
for ((a = 0; a < 4; a++)); do
    for ((b = 0; b < 4; b++)); do
        for ((c = 0; c < 5; c++)); do
            
            # launch our job
            sbatch pso_main_runscript.sh $a $b $c
            sleep 0.5
            
        done        
    done        
done                




