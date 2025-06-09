#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=4             # Total # of nodes 
#SBATCH --ntasks-per-node=32  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=4     # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=03:00:00       # Total run time limit (hh:mm:ss)
#SBATCH -J vort-240612        # Job name
#SBATCH -o vort9.o%j          # Name of stdout output file
#SBATCH -e vort9.e%j          # Name of stderr error file
#SBATCH -p wholenode          # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=all       # Send email at begin and end of job

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
# # srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/evp_three_k0_exec
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec
# 2. add perturbation to the base flow
printf '%s\n' "new_perturb" "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
printf '%s\n' "cor_perturb" "F" "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
# #echo "new_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
# # 3. run initial value code
# echo "new_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/vort9_exec
# echo "new_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/vort_test_exec
echo "cor_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/vort9_exec

# 4. finalize
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
