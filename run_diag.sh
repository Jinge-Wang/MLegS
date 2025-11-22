#!/bin/bash
# SBATCH script to run bsnsq_test on the cluster

#SBATCH -A phy220056          # Allocation name
#SBATCH --nodes=1             # Total # of nodes
#SBATCH --ntasks-per-node=32  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:00
#SBATCH -J test_diag
#SBATCH -o test_diag.o%j
#SBATCH -e test_diag.e%j
#SBATCH -p shared
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

# MAIN SIMULATION EXECUTION:
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_gauss_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_test_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/uz_test_exec
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_critical_exec

# POST-PROCESSING:
# srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/postproc_exec

# PAST TESTS:
# (PASSED) srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/divxm_test_exec
# (PASSED) srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/pe_test_exec
srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/test_shear_exec