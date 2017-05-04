#!/bin/bash

###############################################################################
# predefine parameter
run="run.conf"
NCPU=1
DEBUG_COMMAND_SERIAL="gdb"
DEBUG_COMMAND_PARALLEL_PRE=""
DEBUG_COMMAND_POST_MPIEXEC="xterm -e gdb"

#DEBUG_COMMAND_SERIAL="ddt --connect"
#DEBUG_COMMAND_PARALLEL_PRE=$DEBUG_COMMAND_SERIAL
#DEBUG_COMMAND_POST_MPIEXEC=""

###############################################################################


###############################################################################
# search argument list
for arg in "$@"; do
    par=${arg%%=*}
    val=${arg##*=}
    if [ "$par" == "run" ]; then
        run=$val
    elif [ "$par" == "ncpu" ]; then
        NCPU=$val
    elif [ "$arg" == "-debug" ]; then
        FLAG_DEBUG=1
    elif [ "$arg" == "-no_debugger" ]; then
        FLAG_NO_DEBUGGER=1
    else
        echo "error: unkown parameter " $arg
        exit -1
    fi
done
###############################################################################


###############################################################################
# check for run configuration
if [ ! -f "$run" ]; then
	echo error: missing file "$run"; echo
	exit -1
fi
cp $run run_input
###############################################################################


###############################################################################
# run FLARE

# Absolute path to this script
SCRIPT=$(readlink -f "$0")

# Absolute path this script is in
FLARE_PATH=$(dirname "$SCRIPT")


if [ "$FLAG_DEBUG" == "" ]; then
	if [ "$NCPU" == 1 ]; then
		$FLARE_PATH/flare_bin
	else
		mpiexec -n $NCPU $FLARE_PATH/flare_bin
	fi
else # for debugging only
        if [ "$FLAG_NO_DEBUGGER" == 1 ]; then
	    DEBUG_COMMAND_SERIAL=""
	    DEBUG_COMMAND_PARALLEL_PRE=""
	    DEBUG_COMMAND_POST_MPIEXEC=""
	fi
	if [ "$NCPU" == 1 ]; then
	    $DEBUG_COMMAND_SERIAL $FLARE_PATH/flare_bin_debug
	else
	    $DEBUG_COMMAND_PARALLEL_PRE mpiexec -n $NCPU $DEBUG_COMMAND_POST_MPIEXEC $FLARE_PATH/flare_bin_debug
	fi
fi
###############################################################################


###############################################################################
# cleanup
rm -rf input
rm -rf run_input
###############################################################################
