#!/bin/bash
# SBATCH script to run bsnsq_test on the cluster

#SBATCH -A phy220056          # Allocation name
#SBATCH --nodes=8             # Total # of nodes
#SBATCH --ntasks-per-node=16  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -J bsnsq_test
#SBATCH -o bsnsq_test.o%j
#SBATCH -e bsnsq_test.e%j
#SBATCH -p wholenode
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

# REMOVE ALL OUTPUT FILES
rm -f ./output/*

# --- Automatic Configuration for run_post.sh ---
nodes=$(grep '^#SBATCH --nodes=' $0 | awk -F= '{print $2}' | awk '{print $1}')
ntasks_per_node=$(grep '^#SBATCH --ntasks-per-node=' $0 | awk -F= '{print $2}' | awk '{print $1}')
total_mpi_procs=$((nodes * ntasks_per_node))
if [ "$total_mpi_procs" -gt 128 ]; then
    post_nodes=$(((total_mpi_procs + 127) / 128))
    post_partition="wholenode"
else
    post_nodes=1
    post_partition="shared"
fi
sed -i "s/^#SBATCH --nodes=.*/#SBATCH --nodes=$post_nodes/" run_post.sh
sed -i "s/^#SBATCH --ntasks-per-node=.*/#SBATCH --ntasks-per-node=$((total_mpi_procs / post_nodes))/" run_post.sh
sed -i "s/^#SBATCH -p.*/#SBATCH -p $post_partition/" run_post.sh
# --- End of Automatic Configuration ---

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