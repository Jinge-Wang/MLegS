#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=4   # Total # of MPI tasks per node
#SBATCH --cpus-per-task=4     # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=00:02:00       # Total run time limit (hh:mm:ss)
#SBATCH -J evp_three          # Job name
#SBATCH -o evp_three.o%j      # Name of stdout output file
#SBATCH -e evp_three.e%j      # Name of stderr error file
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
# echo "$SLURM_NTASKS tasks - $SLURM_CPUS_PER_TASK cores per task"

# MPI codes
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/evp_three_centre_exec

# Finalize
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
