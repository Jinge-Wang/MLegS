#!/bin/bash
# Set thread count (default value is 1)...
export OMP_NUM_THREADS=2
ulimit -s unlimited
echo "MPI proc num = "
read NODENUM

# Launch OpenMP code...

# Start Time
start_time="$(date -u +%s.%N)"

# vort2:
# echo "%read.input" | ./bin/copy_exec
# # time echo "%read.input" | ./bin/vort2_exec   # Do not use ibrun or any other MPI launcher
# time ./bin/vort2_exec

# vort:
# make swipe
# echo "%read.input" | ./bin/init_exec

mpirun -np $NODENUM ./bin/init_exec
mpirun -np $NODENUM ./bin/vort_exec

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"
# ---------------------------------------------------

