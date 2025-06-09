#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=8   # Total # of MPI tasks per node
#SBATCH --cpus-per-task=16    # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=02:00:00       # Total run time limit (hh:mm:ss)
#SBATCH -J vort-230130        # Job name
#SBATCH -o vort9.o%j          # Name of stdout output file
#SBATCH -e vort9.e%j          # Name of stderr error file
#SBATCH -p shared             # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=all       # Send email at begin and end of job

# Manage processing environment, load compilers and applications.
module load intel
module load impi
module list

# Set thread count (default value is 1).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Launch MPI code
start_time="$(date -u +%s.%N)"
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec

# # case: cos(mt)cos(kz)
# # Load per_perturb_1 and Add Q-vortex, save to per_perturb_1
# printf '%s\n' "per_perturb_1" "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
# # Load per_perturb_2 and Add to per_perturb_1, save to per_perturb_2
# printf '%s\n' "per_perturb_2" "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_per_exec
# # Run nonlinear simulation
# echo "per_perturb_2" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/vort9_exec

# case: phase adjustment of the degen pair
# Load per_perturb and perform phase adjustments
printf '%s\n' "per_perturb" "120" "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_pha_exec
# Run nonlinear simulation
echo "per_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/vort9_exec

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"

echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
