#!/bin/bash
# FILENAME:  AWV-job_submission

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=2             # Total # of nodes 
#SBATCH --ntasks-per-node=64  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=1     # cpu-cores per task (default value is 1, >1 for multi-threaded tasks)
#SBATCH --time=00:10:00       # Total run time limit (hh:mm:ss)
#SBATCH -J test-240729        # Job name
#SBATCH -o test.o%j           # Name of stdout output file
#SBATCH -e test.e%j           # Name of stderr error file
#SBATCH -p wholenode          # Queue (partition) name
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end       # Send email at begin and end of job

# Manage processing environment, load compilers and applications.
module load intel
module load impi
module list

# Set thread count (default value is 1).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# echo "$SLURM_NTASKS tasks - $SLURM_CPUS_PER_TASK cores per task"

# MPI codes
start_time="$(date -u +%s.%N)"

# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/test4_exec

for PROC_NUM in 4 8 16 32 64 128
do
    echo "$PROC_NUM MPI tasks - $SLURM_CPUS_PER_TASK cores per task"

    srun --mpi=pmi2 -n $PROC_NUM ./bin/test_exchange_exec

    # srun --mpi=pmi2 -n $PROC_NUM ./bin/init_exec
    # printf '%s\n' "new_perturb" "T" | srun --mpi=pmi2 -n $PROC_NUM ./bin/addperturb_non_exec
    # echo "new_perturb" | srun --mpi=pmi2 -n $PROC_NUM ./bin/test2_exec

    # for file in ./test2X/*_GLB_FFF.dat; do mv "$file" "$file.$PROC_NUM"; done
    # for file in ./test2X/*_GLB_PFF.dat; do mv "$file" "$file.$PROC_NUM"; done
    # for file in ./test2X/*_GLB_PPP.dat; do mv "$file" "$file.$PROC_NUM"; done

    # make swipe
done

# echo "$SLURM_NTASKS MPI tasks - $SLURM_CPUS_PER_TASK cores per task"
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/evp_three_k0_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/evp_three_k0_parallel_exec

# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/init_exec
# printf '%s\n' "new_perturb" "T" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/addperturb_non_exec
# echo "new_perturb" | srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/test2_exec

# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/test3_exec

# 4. finalize
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"
echo "GNU Parallel Job Complete"
echo " "
