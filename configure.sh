#!/bin/bash


################################################################################
# predefine parameter

bindir=$HOME/local/bin
libdir=$HOME/local/lib
pythondir=$libdir/python
datadir=$HOME/Database/Magnetic_Configuration

EMC3_dir=""

EXTERNAL_DIR=external
################################################################################


################################################################################
# scanning arguments
for opt in "$@"; do
    par=${opt%%=*}
    val=${opt##*=}
    if [ "$par" == "--help" ]; then
        echo "Configuration script for FLARE"
        echo ""
        echo "Optional arguments are:"
        echo "  --datadir=DIR           data directory for magnetic configurations"
	echo ""
	echo "  --bindir=DIR            directory for installing FLARE executables"
	echo ""
	echo "  --libdir=DIR            directory for the FLARE library"
	echo ""
	echo "  --pythondir=DIR         directory for python modules"
	echo ""
	echo "  --emc3_dir=DIR          set EMC3 source directory"
	echo ""
	echo "  --fusion_io_dir=DIR     set directory of Fusion-IO installation"
	echo ""
	echo "  --with-fgsl             compile with FGSL support"
	echo ""

	exit
    elif [ "$par" == "--datadir" ]; then
        datadir=$val
    elif [ "$par" == "--bindir" ]; then
        bindir=$val
    elif [ "$par" == "--libdir" ]; then
        libdir=$val
    elif [ "$par" == "--pythondir" ]; then
        pythondir=$val
    elif [ "$par" == "--emc3_dir" ]; then
        emc3_dir=$val
    elif [ "$par" == "--fusion_io_dir" ]; then
        fusion_io_dir=$val
    elif [ "$par" == "--with_fgsl" ]; then
        FGSL_SUPPORT=1
    else
        echo "error: unkown parameter " $par
    fi
done
echo
################################################################################


################################################################################
# generate configuration file
rm -rf include.mk
touch include.mk
rm -rf config.h
touch config.h

LOG_FILE=configure.log
echo "command: $0 $@" > $LOG_FILE
echo "timestamp: $(date)" >> $LOG_FILE

# Main program
echo "# Main program"             >> include.mk
echo "PROGRAM        = flare"     >> include.mk
echo "RUN_FLARE      = run_flare.sh" >> include.mk
echo "FLARELIB       = libFLARE.so" >> include.mk
echo "FLAREUI        = flare_ui"  >> include.mk
echo "" >> include.mk


# directory for links to executables
echo "# set up directories for executables and data" >> include.mk
echo "BINDIR         = $bindir" >> include.mk
echo "LIBDIR         = $libdir" >> include.mk
echo "PYTHONDIR      = $pythondir" >> include.mk
echo "Binary directory is $bindir" | tee -a $LOG_FILE
echo "Library directory is $libdir" | tee -a $LOG_FILE
echo "Python module directory is $pythondir" | tee -a $LOG_FILE


# data directory
echo "#define DATABASE_DIR '$datadir'" >> config.h
echo "DATADIR        = $datadir" >> include.mk
echo "Data directory is $datadir" | tee -a $LOG_FILE
echo "" >> include.mk
echo
# ------------------------------------------------------------------------------


# MPI support / Fortran compiler -----------------------------------------------
echo "# Fortran compiler and options"			>> include.mk
if type mpif90 >/dev/null 2>/dev/null; then
	echo "Compiling with MPI support" | tee -a $LOG_FILE
	echo "COMPILER       = mpif90 -DMPI" >> include.mk
elif type ifort >/dev/null 2>/dev/null; then
	echo "Using Intel Fortran compiler" | tee -a $LOG_FILE
	echo "COMPILER       = ifort" >> include.mk
elif type gfortran >/dev/null 2>/dev/null; then
	echo "Using GNU Fortran compiler" | tee -a $LOG_FILE
	echo "COMPILER       = gfortran" >> include.mk
fi
echo "OPT            = -O2 -fconvert=big-endian" >> include.mk
echo "OPT_DEBUG      = -g  -fconvert=big-endian" >> include.mk
echo "" >> include.mk
# ------------------------------------------------------------------------------


# GNU Scientific Library -------------------------------------------------------
if [ "$FGSL_SUPPORT" == "" ]; then
	NOTE="Compiling without FGSL support"
	echo $NOTE | tee -a $LOG_FILE
else
	NOTE="Compiling with FGSL support"
	echo $NOTE | tee -a $LOG_FILE

	FGSL_CFLAGS=`pkg-config --cflags fgsl`
	FGSL_LIBS=`pkg-config --libs fgsl`
	if [ "$FGSL_LIBS" == "" ]; then
		echo "error: cannot find FGSL libraries"
		exit
	fi

	gsl_version=`pkg-config --modversion gsl`
	gsl_version_major=$(echo $gsl_version | cut -d. -f1)
	gsl_version_minor=$(echo $gsl_version | cut -d. -f2)
	gsl_version_micro=$(echo $gsl_version | cut -d. -f3)

	fgsl_version=`pkg-config --modversion fgsl`
	fgsl_version_major=$(echo $fgsl_version | cut -d. -f1)
	fgsl_version_minor=$(echo $fgsl_version | cut -d. -f2)
	fgsl_version_micro=$(echo $fgsl_version | cut -d. -f3)

	echo "#define FGSL" >> config.h
	echo "#define GSL_VERSION_MAJOR_FORTRAN $gsl_version_major" >> config.h
	echo "#define GSL_VERSION_MINOR_FORTRAN $gsl_version_minor" >> config.h

	echo "# GNU Scientific Library"				>> include.mk
	echo "FGSL_FLAGS     = $FGSL_CFLAGS"		        >> include.mk
	echo "FGSL_LIBS      = $FGSL_LIBS"			>> include.mk
	echo ""							>> include.mk
fi
# ------------------------------------------------------------------------------


# set up local source directories ----------------------------------------------
echo "# set up local source directories" >> include.mk
echo "EXTERNAL_DIR   = $EXTERNAL_DIR" >> include.mk
echo "EMC3_LINK_DIR  = \$(EXTERNAL_DIR)/EMC3" >> include.mk
echo "CORE_DIR       = core" >> include.mk
echo "BFIELD_DIR     = bfield" >> include.mk
echo "GEOMETRY_DIR   = geometry" >> include.mk
echo "GRIDGEN_DIR    = fieldline_grid" >> include.mk
echo "TOOLS_DIR      = tools" >> include.mk
echo "DEVEL_DIR      = development" >> include.mk
echo "ADDONS_DIR     = addons" >> include.mk
echo "" >> include.mk
# ------------------------------------------------------------------------------


# checking for coupling to EMC3 ------------------------------------------------
#if [ "$emc3_dir" == "" ]; then
#	NOTE='Compiling without support for fieldline-grid input (based on EMC3 sources)'
#	echo $NOTE | tee -a $LOG_FILE
#else
#	NOTE='Compiling with support for fieldline-grid input (based on EMC3 sources)'
#	echo $NOTE | tee -a $LOG_FILE
#	echo "# $NOTE" >> include.mk
#	echo "EMC3_FLAG      = -DEMC3" >> include.mk
#        echo "EMC3_SRC_DIR   = $emc3_dir/MAIN" >> include.mk
#	if [ ! -d "$emc3_dir/MAIN" ]; then
#		echo "error: EMC3 source files not found!"
#		exit
#	fi
#	echo "EMC3_OBJ       = \
#               PHYS_CONST.o\
#               GEOMETRY_PL.o\
#               SURFACE_PL.o\
#               MAPPING.o\
#               check.o\
#               ibm_iface.o\
#               random.o\
#               real_to_ft.o\
#               service.o\
#               sf_def_user.o\
#               sf_jump.o" >> include.mk
#	echo "EMC3_OBJ_LONG  = \$(addprefix \$(EMC3_LINK_DIR)/,\$(EMC3_OBJ))" >> include.mk
#	echo "" >> include.mk
#
#
#fi
# ------------------------------------------------------------------------------


# checking for Fusion-IO installation ------------------------------------------
if [ "$fusion_io_dir" == "" ]; then
	NOTE='Compiling without Fusion-IO support'
	echo $NOTE | tee -a $LOG_FILE
else
	NOTE='Compiling with Fusion-IO support'

	# check HDF5 installation
	HDF5="hdf5-openmpi"
	HDF5_CFLAGS=`pkg-config --cflags $HDF5`
	HDF5_LIBS=`pkg-config --libs $HDF5`
	echo "HDF5_FLAGS     = $HDF5_CFLAGS"		        >> include.mk
	echo "HDF5_LIBS      = $HDF5_LIBS"			>> include.mk

	echo $NOTE | tee -a $LOG_FILE
	echo "# $NOTE" >> include.mk
	echo "FIO_FLAGS      = -DFIO -I $fusion_io_dir/include" >> include.mk
	echo "FIO_LIBS       = -L $fusion_io_dir/lib -lfusionio -lm3dc1 -lstdc++" >> include.mk
	echo "" >> include.mk
fi
# ------------------------------------------------------------------------------


# ODE solvers ------------------------------------------------------------------
echo "# additional ODE solvers" >> include.mk
declare -a ode_solver_files=("dlsode.f")
declare -a ode_solver_names=("DLSODE"  )
n_ode=${#ode_solver_files[@]}
ODE_FLAGS=''
ODE_OBJECTS=''

for (( i=0; i<${n_ode}; i++ )); do
	if [ -f src/$EXTERNAL_DIR/${ode_solver_files[$i]} ]; then
		NOTE="Compiling with ${ode_solver_names[$i]}"
		echo $NOTE | tee -a $LOG_FILE
		ODE_FLAGS="$ODE_FLAGS -D${ode_solver_names[$i]}"
		ODE_OBJECTS="$ODE_OBJECTS ${ode_solver_files[$i]%%.f*}.o"
	else
		NOTE="Compiling without ${ode_solver_names[$i]}"
		echo $NOTE | tee -a $LOG_FILE
	fi
done
echo "ODE_FLAGS      = $ODE_FLAGS"			>> include.mk
echo "ODE_OBJECTS    = $ODE_OBJECTS"			>> include.mk
echo "" >> include.mk
# ------------------------------------------------------------------------------


# Flags and libraries
FLAGS='$(FGSL_FLAGS) $(HDF5_FLAGS) $(FIO_FLAGS) $(ODE_FLAGS)'
LIBS='$(FGSL_LIBS)  $(HDF5_LIBS)  $(FIO_LIBS)'
echo "# Flags and libraries"				>> include.mk
echo "FLAGS          = -fPIC $FLAGS"			>> include.mk
echo "LIBS           = $LIBS"				>> include.mk
echo ""							>> include.mk


#
echo 'FC             = $(COMPILER) $(FLAGS) -DFLARE $(OPT)' >> include.mk
echo 'FC_ADDONS      = $(COMPILER) $(FLAGS)         $(OPT)' >> include.mk
echo 'FC_DEBUG       = $(COMPILER) $(FLAGS) -DFLARE $(OPT_DEBUG)' >> include.mk
################################################################################
echo
