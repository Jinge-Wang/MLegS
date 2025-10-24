#!/bin/bash
# SBATCH script to run bsnsq_test on the cluster

#SBATCH -A phy220056          # Allocation name
#SBATCH --nodes=1             # Total # of nodes
#SBATCH --ntasks-per-node=64  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=2
#SBATCH --time=01:30:00
#SBATCH -J bsnsq_test
#SBATCH -o bsnsq_test.o%j
#SBATCH -e bsnsq_test.e%j
#SBATCH -p shared
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

# REMOVE ALL OUTPUT FILES
rm -f ./output/*

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "$SLURM_NTASKS tasks - $SLURM_CPUS_PER_TASK cores per task"

# MAIN SIMULATION EXECUTION:
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_gauss_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_test_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/uz_test_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_critical_exec

# POST-PROCESSING:
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/postproc_exec

# PAST TESTS:
# (PASSED) srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/divxm_test_exec
# (PASSED) srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/pe_test_exec