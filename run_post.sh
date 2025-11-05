#!/bin/bash
# SBATCH script to run postproc on the cluster

#SBATCH -A phy220056          # Allocation name 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00
#SBATCH -J post
#SBATCH -o postproc.o%j
#SBATCH -e postproc.e%j
#SBATCH -p shared
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

# Manage processing environment, load compilers and applications.
module load intel
module load impi
module list

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "$SLURM_NTASKS tasks - $SLURM_CPUS_PER_TASK cores per task"

srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/postproc_exec
