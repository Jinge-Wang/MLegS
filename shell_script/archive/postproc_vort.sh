#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=8   # Total # of MPI tasks per node
#SBATCH --cpus-per-task=16    # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=00:02:00       # Total run time limit (hh:mm:ss)
#SBATCH -J postproc_vort      # Job name
#SBATCH -o postproc_vort.o    # Name of stdout output file
#SBATCH -e postproc_vort.e    # Name of stderr error file
#SBATCH -p debug              # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=none      # Send email at begin and end of job

# Manage processing environment, load compilers and applications.
module load intel
module load impi
module list

# Set thread count (default value is 1).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Launch MPI code
start_time="$(date -u +%s.%N)"
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/postproc_vort_exec
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"

echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
