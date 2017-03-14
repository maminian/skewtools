#!/bin/bash
# Generate a picture profiling a previously run script.
#
# Piping sexiness due to:
#
# https://benjimenez.wordpress.com/2012/02/19/profiling-a-simple-fortran-code-with-gprof/

# Pull the directory this shell script is located; it should also have gprof2dot.py.
# First get the script's path...
SCRIPT=$(readlink -f $0)
# then get the directory containing the script.
SCRIPTPATH=$(dirname "$SCRIPT")

#echo $SCRIPTPATH/gprof2dot.py

# Finally, do the call.
# The options for gprof2dot.py force it to display the entire evaluation graph.
gprof $1 | "$SCRIPTPATH/gprof2dot.py" --node-thres=0.01 --edge-thres=0.01 | dot -Tpng -o profile.png
eog profile.png &
