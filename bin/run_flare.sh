#!/bin/bash

###############################################################################
# predefine parameter
procs=1
run_default="run.conf"

DEBUG_COMMAND_SERIAL="gdb --args"
DEBUG_COMMAND_PARALLEL_PRE=""
DEBUG_COMMAND_POST_MPIEXEC="xterm -e gdb --args"

#DEBUG_COMMAND_SERIAL="ddt --connect"
#DEBUG_COMMAND_PARALLEL_PRE=$DEBUG_COMMAND_SERIAL
#DEBUG_COMMAND_POST_MPIEXEC=""
###############################################################################


###############################################################################
# run control
arg_list="$1"

# provide user defined run control file
if [ "$1" == "run" ]; then
    shift
    if [ $# -eq 0 ]; then
        echo "error: missing argument for run control file!"
        exit -1
    fi

    run_conf="$1"
    if [ ! -f "$run_conf" ]; then
        echo "error: run control file $run_conf does not exist!"
        exit -1
    fi
    arg_list="$arg_list $run_conf"
    shift


# import equilibrium into database
elif [ "$1" == "import" ]; then
    shift
    if [ $# -eq 0 ]; then
        echo "error: filename missing for import equilibrium!"
        exit -1
    fi
    arg_list="$arg_list $1"
    shift


# use default run control file
else
    arg_list="run $run_default"
fi


# browse through remaining argument list
while test $# -gt 0
do
    # parallel execution
    if [ "$1" == "parallel" ]; then
        shift
        if [ ! -z "${1##*[!0-9]*}" ]; then
            procs=$1
        else
            echo "error: number of processes missing for argument 'parallel'!"
            exit -1
        fi


    # debug mode
    elif [ "$1" == "debug" ]; then
        FLAG_DEBUG=1


    # undefined argument
    else
        echo "error: undefined argument $1"
        exit -1
    fi

    shift
done
###############################################################################


###############################################################################
# run FLARE

# Absolute path to this script
SCRIPT=$(readlink -f "$0")

# Absolute path this script is in
FLARE_PATH=$(dirname "$SCRIPT")


if [ "$FLAG_DEBUG" == "" ]; then
	if [ "$procs" == 1 ]; then
		$FLARE_PATH/FLARE $arg_list
	else
		mpiexec -n $procs $FLARE_PATH/FLARE $arg_list
	fi
else # for debugging only
	if [ "$procs" == 1 ]; then
	    $DEBUG_COMMAND_SERIAL $FLARE_PATH/FLARE $arg_list
	else
	    $DEBUG_COMMAND_PARALLEL_PRE mpiexec -n $procs $DEBUG_COMMAND_POST_MPIEXEC $FLARE_PATH/FLARE $arg_list
	fi
fi
###############################################################################
