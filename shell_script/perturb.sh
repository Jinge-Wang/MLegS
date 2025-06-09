#!/bin/bash
# Set thread count (default value is 1)...
export OMP_NUM_THREADS=4

# Launch OpenMP code...
make swipe
echo "%read.input" | ./bin/init_exec
echo "%read.input" | ./bin/copy_exec
echo "%read.input" | ./bin/addperturb_exec
echo "%read.input" | ./bin/copy_exec
 echo "%read.input" | ./bin/vort_exec
#./bin/vort2_exec


# ---------------------------------------------------

