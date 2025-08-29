#!/bin/bash
# SBATCH script to run bsnsq_test on the cluster

#SBATCH -A phy220056          # Allocation name
#SBATCH --nodes=1             # Total # of nodes
#SBATCH --ntasks-per-node=32  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH -J bsnsq_test
#SBATCH -o bsnsq_test.o%j
#SBATCH -e bsnsq_test.e%j
#SBATCH -p shared
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# REMOVE ALL OUTPUT FILES
rm -f ./output/*

echo "Running bsnsq_test with $SLURM_NTASKS MPI tasks"

srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_test_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/uz_test_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_critical_exec

# PAST TESTS:
# (PASSED) srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/divxm_test_exec
# (PASSED) srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/pe_test_exec