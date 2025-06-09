#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=2             # Total # of nodes 
#SBATCH --ntasks-per-node=64  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=2     # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=00:10:00       # Total run time limit (hh:mm:ss)
#SBATCH -J post-241006        # Job name
#SBATCH -o postproc.o%j       # Name of stdout output file
#SBATCH -e postproc.e%j       # Name of stderr error file
#SBATCH -p debug             # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end       # Send email at begin and end of job

# Manage processing environment, load compilers and applications.
module load intel
module load impi
module list

# Set thread count (default value is 1).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "$SLURM_NTASKS tasks - $SLURM_CPUS_PER_TASK cores per task"

# MPI codes
# 1. initialize the base flow
start_time="$(date -u +%s.%N)"
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/postproc_mpi_exec

# 2. finalize
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
