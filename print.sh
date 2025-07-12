#!/bin/bash

# Script to run evp_print_exec for various combinations of m and k
# Usage: ./print.sh

# Check if the executable exists
if [ ! -f "./bin/evp_print_exec" ]; then
    echo "Error: ./bin/evp_print_exec not found"
    exit 1
fi

# Set OMP
export OMP_NUM_THREADS=4

# Loop through m values (0, 1, 2)
for m in 0 1 2; do
    # Loop through k values (-2, -1, 0, 1, 2)
    for k in -2 -1 0 1 2; do
        echo "Running with m=$m, k=$k"
        mpirun -np 4 ./bin/evp_print_exec $m $(printf "%.1f" $k)
        echo "------------------------"
    done
done

echo "All runs completed."