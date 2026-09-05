#!/bin/bash
# Parameter scanning script for evp_print_exec
# Supports both local execution and SLURM supercomputer submission

# SLURM Configuration (used when running on supercomputers)
SLURM_ACCOUNT="phy220056"
SLURM_NODES=1
SLURM_NTASKS_PER_NODE=8
SLURM_CPUS_PER_TASK=4
SLURM_TIME="00:15:00"
SLURM_PARTITION="shared"
SLURM_EMAIL="jinge@berkeley.edu"

# Default parameter values
BV=0.0 #$(echo "1/0.9" | bc -l)  # Brunt-Väisälä frequency using arithmetic
OMEGA=0.0     # Angular velocity
NRCHOP=350    # Radial resolution
M_START=0     # Start value for m
M_END=3       # End value for m
K_START=0.0  # Start value for k
K_END=10.0     # End value for k
K_STEP=0.025    # Step size for k

# System detection
FORCE_LOCAL=false
FORCE_SLURM=false

# Function to detect if we're on a SLURM system
detect_slurm_system() {
    if command -v sbatch &> /dev/null && command -v srun &> /dev/null; then
        return 0  # SLURM system detected
    else
        return 1  # Not a SLURM system
    fi
}

# Function to print usage
print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  --bv=VALUE        Brunt-Väisälä frequency (default: $BV)"
    echo "  --omega=VALUE     Angular velocity (default: $OMEGA)"
    echo "  --nr=VALUE        Radial resolution (default: $NRCHOP)"
    echo "  --mstart=VALUE    Start value for m (default: $M_START)"
    echo "  --mend=VALUE      End value for m (default: $M_END)"
    echo "  --kstart=VALUE    Start value for k (default: $K_START)"
    echo "  --kend=VALUE      End value for k (default: $K_END)"
    echo "  --kstep=VALUE     Step size for k (default: $K_STEP)"
    echo "  --local           Force local execution (ignore SLURM detection)"
    echo "  --slurm           Force SLURM execution"
    echo "  --help            Show this help message"
}

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
    --local)
      FORCE_LOCAL=true
      shift
      ;;
    --slurm)
      FORCE_SLURM=true
      shift
      ;;
    --help)
      print_usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Determine execution mode
if [ "$FORCE_LOCAL" = true ]; then
    USE_SLURM=false
    echo "Forced local execution mode"
elif [ "$FORCE_SLURM" = true ]; then
    USE_SLURM=true
    echo "Forced SLURM execution mode"
elif detect_slurm_system; then
    USE_SLURM=true
    echo "SLURM system detected - using SLURM execution mode"
else
    USE_SLURM=false
    echo "Local system detected - using local execution mode"
fi

# Determine executable name based on parameters
# If both bv and omega are 0, use the original executable
if [ $(echo "$BV == 0.0" | bc -l) -eq 1 ] && [ $(echo "$OMEGA == 0.0" | bc -l) -eq 1 ]; then
    EXECUTABLE="./bin/evp_print_org_exec"
    USE_BV_OMEGA=false
    echo "Both bv and omega are 0 - using original executable (no stratification/rotation)"
else
    EXECUTABLE="./bin/evp_print_exec"
    USE_BV_OMEGA=true
    echo "Using Boussinesq executable with bv=$BV, omega=$OMEGA"
fi

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: $EXECUTABLE not found"
    echo "Please make sure the executable is compiled and available"
    exit 1
fi

echo "Using executable: $EXECUTABLE"

# Create output directory if it doesn't exist
mkdir -p data

# # Remove old output files
# if ls data/bsnsq_eig_bv_* 1> /dev/null 2>&1; then
#     echo "Removing old output files in data/ directory"
#     rm data/bsnsq_eig_bv_*
# fi

# Display run parameters
echo "==================== Parameter Sweep ======================"
echo "Execution mode: $([ "$USE_SLURM" = true ] && echo "SLURM" || echo "Local")"
echo "Executable: $EXECUTABLE"
if [ "$USE_BV_OMEGA" = true ]; then
    echo "Brunt-Väisälä (bv): $BV"
    echo "Angular velocity (w): $OMEGA"
else
    echo "No stratification/rotation (bv=0, w=0)"
fi
echo "Radial resolution (nr): $NRCHOP"
echo "m range: $M_START to $M_END"
echo "k range: $K_START to $K_END (step: $K_STEP)"
echo "=========================================================="

# Function to run the job locally
run_local() {
    local m=$1
    local k_formatted=$2
    
    # Set OMP threads for local execution
    export OMP_NUM_THREADS=3
    
    echo "Running locally with mpirun -np 4..."
    if [ "$USE_BV_OMEGA" = true ]; then
        mpirun -np 4 $EXECUTABLE m=$m k=$k_formatted bv=$BV w=$OMEGA nr=$NRCHOP
    else
        mpirun -np 4 $EXECUTABLE m=$m k=$k_formatted nr=$NRCHOP
    fi
    return $?
}

# Function to run the job via SLURM
run_slurm() {
    # Create a single SLURM script that handles the entire parameter sweep
    local script_name="slurm_parameter_sweep.sh"
    
    cat > "$script_name" << EOF
#!/bin/bash
#SBATCH -A $SLURM_ACCOUNT
#SBATCH --nodes=$SLURM_NODES
#SBATCH --ntasks-per-node=$SLURM_NTASKS_PER_NODE
#SBATCH --cpus-per-task=$SLURM_CPUS_PER_TASK
#SBATCH --time=$SLURM_TIME
#SBATCH -J evp_print_sweep
#SBATCH -o print_sweep.o%j
#SBATCH -e print_sweep.e%j
#SBATCH -p $SLURM_PARTITION
#SBATCH --mail-user=$SLURM_EMAIL
#SBATCH --mail-type=end

# Load modules
module load intel
module load impi

# Set thread count
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

echo "Starting parameter sweep on SLURM cluster"
if [ "$USE_BV_OMEGA" = true ]; then
    echo "Parameters: BV=$BV, OMEGA=$OMEGA, NRCHOP=$NRCHOP"
else
    echo "Parameters: NRCHOP=$NRCHOP (no stratification/rotation)"
fi
echo "m range: $M_START to $M_END"
echo "k range: $K_START to $K_END (step: $K_STEP)"

# Counter for completed runs
total_runs=0
skipped_runs=0

# Loop through m values
for m in \$(seq $M_START $M_END); do
    # Loop through k values using bc for floating point arithmetic
    k=$K_START
    while (( \$(echo "\$k <= $K_END" | bc -l) )); do
        # Format k value for display
        k_formatted=\$(printf "%.2f" \$k)
        
        # Skip the case m=0, k=0
        if [ "\$m" -eq "0" ] && [ \$(echo "\$k == 0.0" | bc -l) -eq 1 ]; then
            echo "Skipping m=0, k=0 (trivial case)"
            skipped_runs=\$((skipped_runs + 1))
        else
            echo "==============================================="
            if [ "$USE_BV_OMEGA" = true ]; then
                echo "Running with m=\$m, k=\$k_formatted, bv=$BV, w=$OMEGA, nr=$NRCHOP"
            else
                echo "Running with m=\$m, k=\$k_formatted, nr=$NRCHOP"
            fi
            echo "==============================================="
            
            # Run the eigenvalue calculation
            if [ "$USE_BV_OMEGA" = true ]; then
                srun --mpi=pmi2 -n \$SLURM_NTASKS $EXECUTABLE m=\$m k=\$k_formatted bv=$BV w=$OMEGA nr=$NRCHOP
            else
                srun --mpi=pmi2 -n \$SLURM_NTASKS $EXECUTABLE m=\$m k=\$k_formatted nr=$NRCHOP
            fi
            
            # Check if the run was successful
            if [ \$? -eq 0 ]; then
                echo "Run completed successfully"
                total_runs=\$((total_runs + 1))
            else
                echo "Run failed!"
            fi
        fi
        
        # Increment k by step size
        k=\$(echo "\$k + $K_STEP" | bc -l)
    done
done

echo "=========================================================="
echo "Parameter sweep completed on SLURM:"
echo "- Total successful runs: \$total_runs"
echo "- Skipped combinations: \$skipped_runs"
echo "- Output files are in the data/ directory"
EOF

    echo "Submitting single SLURM job for parameter sweep: $script_name"
    sbatch "$script_name"
    local submit_status=$?
    
    # Keep the script file for reference
    echo "SLURM script saved as: $script_name"
    
    return $submit_status
}

# Counter for completed runs
total_runs=0
skipped_runs=0

# Execute based on system type
if [ "$USE_SLURM" = true ]; then
    # For SLURM systems, submit one job that handles the entire parameter sweep
    echo "Submitting SLURM job for parameter sweep..."
    run_slurm
    if [ $? -eq 0 ]; then
        echo "SLURM job submitted successfully"
        echo "Check job status with: squeue -u \$USER"
        echo "Monitor output with: tail -f print_sweep.o<job_id>"
    else
        echo "Failed to submit SLURM job!"
        exit 1
    fi
else
    # For local systems, run the parameter sweep directly
    echo "Running parameter sweep locally..."
    
    # Loop through m values
    for m in $(seq $M_START $M_END); do
        # Loop through k values
        k=$K_START
        while (( $(echo "$k <= $K_END" | bc -l) )); do
            # Format k value for display
            k_formatted=$(printf "%.3f" $k)
            
            # Skip the case m=0, k=0
            if [ "$m" -eq "0" ] && [ $(echo "$k == 0.0" | bc -l) -eq 1 ]; then
                echo "Skipping m=0, k=0 (trivial case)"
                skipped_runs=$((skipped_runs + 1))
            else
                echo "==============================================="
                if [ "$USE_BV_OMEGA" = true ]; then
                    echo "Running with m=$m, k=$k_formatted, bv=$BV, w=$OMEGA, nr=$NRCHOP"
                else
                    echo "Running with m=$m, k=$k_formatted, nr=$NRCHOP"
                fi
                echo "==============================================="
                
                # Run the eigenvalue calculation locally
                run_local $m $k_formatted
                run_status=$?
                
                # Check if the run was successful
                if [ $run_status -eq 0 ]; then
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
fi