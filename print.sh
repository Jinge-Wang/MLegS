#!/bin/bash
# filepath: /Users/jinge/Projects/MLegS-dev/MLegS/print.sh
# Parameter scanning script for evp_print_exec

# Default parameter values
BV=0.0 #$(echo "1/0.9" | bc -l)  # Brunt-Väisälä frequency using arithmetic
OMEGA=-0.6     # Angular velocity
NRCHOP=350    # Radial resolution
M_START=0     # Start value for m
M_END=0       # End value for m
K_START=0.0  # Start value for k
K_END=50.0     # End value for k
K_STEP=0.5    # Step size for k

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --bv=*)
      BV="${1#*=}"
      shift
      ;;
    --omega=*)
      OMEGA="${1#*=}"
      shift
      ;;
    --nr=*)
      NRCHOP="${1#*=}"
      shift
      ;;
    --mstart=*)
      M_START="${1#*=}"
      shift
      ;;
    --mend=*)
      M_END="${1#*=}"
      shift
      ;;
    --kstart=*)
      K_START="${1#*=}"
      shift
      ;;
    --kend=*)
      K_END="${1#*=}"
      shift
      ;;
    --kstep=*)
      K_STEP="${1#*=}"
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check if the executable exists
if [ ! -f "./bin/evp_print_exec" ]; then
    echo "Error: ./bin/evp_print_exec not found"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p data

# # Remove old output files
# if ls data/bsnsq_eig_bv_* 1> /dev/null 2>&1; then
#     echo "Removing old output files in data/ directory"
#     rm data/bsnsq_eig_bv_*
# fi

# Display run parameters
echo "==================== Parameter Sweep ======================"
echo "Brunt-Väisälä (bv): $BV"
echo "Angular velocity (w): $OMEGA"
echo "Radial resolution (nr): $NRCHOP"
echo "m range: $M_START to $M_END"
echo "k range: $K_START to $K_END (step: $K_STEP)"
echo "=========================================================="

# Set OMP threads
export OMP_NUM_THREADS=3

# Counter for completed runs
total_runs=0
skipped_runs=0

# Loop through m values
for m in $(seq $M_START $M_END); do
    # Loop through k values
    k=$K_START
    while (( $(echo "$k <= $K_END" | bc -l) )); do
        # Format k value for display
        k_formatted=$(printf "%.2f" $k)
        
        # Skip the case m=0, k=0
        if [ "$m" -eq "0" ] && [ $(echo "$k == 0.0" | bc -l) -eq 1 ]; then
            echo "Skipping m=0, k=0 (trivial case)"
            skipped_runs=$((skipped_runs + 1))
        else
            echo "==============================================="
            echo "Running with m=$m, k=$k_formatted, bv=$BV, w=$OMEGA, nr=$NRCHOP"
            echo "==============================================="
            
            # Run the eigenvalue calculation with named parameters
            mpirun -np 4 ./bin/evp_print_exec m=$m k=$k_formatted bv=$BV w=$OMEGA nr=$NRCHOP
            
            # Check if the run was successful
            if [ $? -eq 0 ]; then
                echo "Run completed successfully"
                total_runs=$((total_runs + 1))
            else
                echo "Run failed!"
            fi
        fi
        
        # Increment k by step size
        k=$(echo "$k + $K_STEP" | bc -l)
    done
done

echo "=========================================================="
echo "Parameter sweep completed:"
echo "- Total successful runs: $total_runs"
echo "- Skipped combinations: $skipped_runs"
echo "- Output files are in the data/ directory"