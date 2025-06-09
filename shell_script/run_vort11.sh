#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=8   # Total # of MPI tasks per node
#SBATCH --cpus-per-task=16    # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=02:00:00       # Total run time limit (hh:mm:ss)
#SBATCH -J vort-230330        # Job name
#SBATCH -o vort11.o%j         # Name of stdout output file
#SBATCH -e vort11.e%j         # Name of stderr error file
#SBATCH -p shared             # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=all       # Send email at begin and end of job

# Manage processing environment, load compilers and applications.
module load intel
module load impi
module list

# Set thread count (default value is 1).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Perturbation file name
filename="per_perturb"

# Launch MPI code
start_time="$(date -u +%s.%N)"
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec
echo "$filename" | ./bin/addperturb_split_exec
printf '%s\n' "$filename"0 "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
printf '%s\n' "$filename"1 "F" "F"| srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
#echo "new_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
echo "$filename" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/vort11_exec
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"

echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
