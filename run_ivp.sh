#!/bin/bash
# IVP (Initial Value Problem) job submission script
# Runs initialization (init_exec) followed by Boussinesq IVP simulation (bsnsq_ivp_exec)

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=16  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=2     # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=00:05:00       # Total run time limit (hh:mm:ss) - longer for IVP
#SBATCH -J bsnsq_ivp          # Job name
#SBATCH -o ivp.o%j            # Name of stdout output file
#SBATCH -e ivp.e%j            # Name of stderr error file
#SBATCH -p shared             # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end       # Send email at begin and end of job

# Manage processing environment, load compilers and applications.
start_time="$(date -u +%s.%N)"
module load intel
module load impi
module list

# Set thread count (default value is 1).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Running with $SLURM_NTASKS MPI tasks - $SLURM_CPUS_PER_TASK cores per task"

# Print job information
echo "=========================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Start Time: $(date)"
echo "Working Directory: $(pwd)"
echo "=========================================================="

# Check if required executables exist
if [ ! -f "./bin/init_exec" ]; then
    echo "Error: ./bin/init_exec not found"
    exit 1
fi

if [ ! -f "./bin/bsnsq_ivp_exec" ]; then
    echo "Error: ./bin/bsnsq_ivp_exec not found"
    exit 1
fi

# Check if input file exists
if [ ! -f "read.input" ]; then
    echo "Error: read.input file not found"
    exit 1
fi

echo "All required files found. Starting IVP simulation..."

# Step 1: Run initialization
echo "=========================================================="
echo "Step 1: Running initialization (init_exec)"
echo "Start time: $(date)"
echo "=========================================================="

srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec

init_status=$?
if [ $init_status -ne 0 ]; then
    echo "Error: Initialization failed with exit code $init_status"
    exit $init_status
fi

echo "Initialization completed successfully"

# Step 2: Run Boussinesq IVP simulation
echo "=========================================================="
echo "Step 2: Running Boussinesq IVP simulation (bsnsq_ivp_exec)"
echo "Start time: $(date)"
echo "=========================================================="

# Provide empty input to use default filename when prompted
echo "" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_ivp_exec

ivp_status=$?
if [ $ivp_status -ne 0 ]; then
    echo "Error: IVP simulation failed with exit code $ivp_status"
    exit $ivp_status
fi

echo "IVP simulation completed successfully"

# Finalize and report timing
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"

echo "=========================================================="
echo "Job completed successfully!"
echo "End time: $(date)"
echo "Total elapsed time: ${elapsed} seconds"
echo "Check output files in the data/ or output/ directory"
echo "=========================================================="
