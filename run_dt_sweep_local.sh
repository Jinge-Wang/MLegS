#!/bin/bash
# Local script to run bsnsq_gauss_exec over a sweep of time step sizes
# Uses mpirun for execution (no SLURM environment required)

# ======================================================================
# MPI configuration – adjust these for your local machine
# ======================================================================
MPI_PROCS=4
OMP_THREADS=1

# ======================================================================
# Time step sizes to sweep over
# ======================================================================
DT_VALUES=(0.1 0.05 0.025 0.01)

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

export OMP_NUM_THREADS=$OMP_THREADS
echo "Using $MPI_PROCS MPI processes, $OMP_THREADS OpenMP threads per process."

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
    mpirun -np $MPI_PROCS ./bin/bsnsq_gauss_exec

    echo "  Run with DT = ${dt} complete."

done

echo "======================================================"
echo "DT sweep complete. All runs finished."
echo "======================================================"
