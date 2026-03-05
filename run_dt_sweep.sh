#!/bin/bash
# SBATCH script to run bsnsq_gauss_exec over a sweep of time step sizes

#SBATCH -A phy220056          # Allocation name
#SBATCH --nodes=8             # Total # of nodes
#SBATCH --ntasks-per-node=16  # Total # of MPI tasks per node
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -J bsnsq_dt_sweep
#SBATCH -o bsnsq_dt_sweep.o%j
#SBATCH -e bsnsq_dt_sweep.e%j
#SBATCH -p wholenode
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "$SLURM_NTASKS tasks - $SLURM_CPUS_PER_TASK cores per task"

# ======================================================================
# Time step sizes to sweep over
# ======================================================================
DT_VALUES=(0.1 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025)

# ======================================================================
# Helper: modify the DT line in read.input
#   Finds the keyword "DT", skips comment lines (starting with #),
#   then replaces the first number (dt) on the data line while
#   preserving the second number (total simulation time).
# ======================================================================
set_dt() {
    local new_dt=$1
    awk -v OFS='  ' -v new_dt="$new_dt" '
    BEGIN { in_dt=0 }
    /^DT$/ { in_dt=1; print; next }
    in_dt && /^[[:space:]]*#/ { print; next }
    in_dt { $1 = new_dt; in_dt = 0 }
    { print }
    ' read.input > read.input.tmp && mv read.input.tmp read.input
}

# ======================================================================
# REMOVE ALL OUTPUT FILES before the sweep begins (once only)
# ======================================================================
rm -f ./output/*
echo "Output directory cleared."

# ======================================================================
# DT SWEEP LOOP
# ======================================================================
for dt in "${DT_VALUES[@]}"; do

    echo "======================================================"
    echo "  Running bsnsq_gauss_exec with DT = ${dt}"
    echo "======================================================"

    # Modify read.input with the new time step
    set_dt "$dt"
    echo "  read.input updated: DT = ${dt}"

    # Run the simulation
    srun --mpi=pmi2 -n $SLURM_NTASKS ./bin/bsnsq_gauss_exec

    echo "  Run with DT = ${dt} complete."

done

echo "======================================================"
echo "DT sweep complete. All runs finished."
echo "======================================================"
