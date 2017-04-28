#!/bin/bash

###############################################################################
# predefine parameter
procs=1
run_default="run.conf"
###############################################################################


###############################################################################
# run control
run_conf="$2"
if [ "$1" == "run" ]; then
    shift
    if [ $# -eq 0 ]; then
        echo "error: missing argument for run control file!"
        exit -1
    fi
    shift

elif [ "$1" == "import" ]; then
    shift
    if [ $# -eq 0 ]; then
        echo "error: filename missing for import equilibrium!"
        exit -1
    fi
    shift

else
    run_conf=$run_default
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
		$FLARE_PATH/flare_bin $run_conf
	else
		mpiexec -n $procs $FLARE_PATH/flare_bin $run_conf
	fi
else # for debugging only
	if [ "$procs" == 1 ]; then
		gdb $FLARE_PATH/flare_bin_debug $run_conf
	else
		mpiexec -n $procs xterm -e gdb $FLARE_PATH/flare_bin_debug $run_conf
	fi
fi
###############################################################################
