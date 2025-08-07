#!/bin/bash
# SBATCH script to run test_prodct on the cluster

#SBATCH -A phy220056          # Allocation name
#SBATCH --nodes=1             # Total # of nodes
#SBATCH --ntasks-per-node=4   # Total # of MPI tasks per node
#SBATCH --cpus-per-task=4
#SBATCH --time=00:02:00
#SBATCH -J test_prodct
#SBATCH -o test_prodct.o%j
#SBATCH -e test_prodct.e%j
#SBATCH -p shared
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running test_prodct with $SLURM_NTASKS MPI tasks"

srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/test_prodct_exec