#!/bin/bash


################################################################################
# predefine parameter
base_dir="Database/Magnetic_Configuration"

EMC3_dir=""

fusion_io_dir=$M3DC1_INSTALL_DIR
if [ "$M3DC1_ARCH" == "" ]; then
    fusion_io_arch=`uname`
else
    fusion_io_arch=$M3DC1_ARCH
fi


BIN_DIR=$HOME/bin

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
        echo "  --base_dir=DIR          base directory for magnetic configuration"
	echo "	                        database (relative to home directory)"
	echo ""
	echo "  --bin_dir=DIR           directory with links to executables"
	echo ""
	echo "  --emc3_dir=DIR          set EMC3 source directory"
	echo ""
	echo "  --fusion_io_dir=DIR     set directory of M3D-C1 / fusion_io installation"
	echo "  --fusion_io_arch=ARCH   set architecture string used for M3D-C1 (if not equal to `uname`)"
	echo ""

	exit
    elif [ "$par" == "--base_dir" ]; then
        base_dir=$val
    elif [ "$par" == "--bin_dir" ]; then
        BIN_DIR=$val
    elif [ "$par" == "--emc3_dir" ]; then
        emc3_dir=$val
    elif [ "$par" == "--fusion_io_dir" ]; then
        fusion_io_dir=$val
    elif [ "$par" == "--fusion_io_arch" ]; then
        fusion_io_arch=$val
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

LOG_FILE=configure.log
echo "command: $0 $@" > $LOG_FILE
echo "timestamp: $(date)" >> $LOG_FILE

# Main program
echo "# Main program"             >> include.mk
echo "PROGRAM        = flare_bin" >> include.mk
echo "PROGRAM_DEBUG  = flare_bin_debug" >> include.mk
echo "" >> include.mk


# directory for links to executables
echo "# setting directory for links to executables" >> include.mk
echo "BIN_DIR        = $BIN_DIR" >> include.mk
echo "" >> include.mk
echo "Binary directory is $BIN_DIR" | tee -a $LOG_FILE
# ------------------------------------------------------------------------------


# MPI support / Fortran compiler -----------------------------------------------
echo "# Fortran compiler and options"			>> include.mk
if type mpif90 >/dev/null 2>/dev/null; then
	echo "Compiling with MPI support" | tee -a $LOG_FILE
	echo "COMPILER       = mpif90 -DMPI" >> include.mk
elif type ifort >/dev/null 2>/dev/null; then
	echo "Using Inter Fortran compiler" | tee -a $LOG_FILE
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
FGSL_CFLAGS=`pkg-config --cflags fgsl`
FGSL_LIBS=`pkg-config --libs fgsl`
if [ "$FGSL_LIBS" == "" ]; then
	echo "Warning: FGSL not found! Some tools may not work properly!"
fi
echo "# GNU Scientific Library"				>> include.mk
echo "FGSL_CFLAGS    = $FGSL_CFLAGS"			>> include.mk
echo "FGSL_LIBS      = $FGSL_LIBS"			>> include.mk
echo ""							>> include.mk
# ------------------------------------------------------------------------------


# setting local source directories ---------------------------------------------
echo "# setting local source directories" >> include.mk
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


# checking for M3D-C1/fusion_io installation -----------------------------------
if [ "$fusion_io_dir" == "" ]; then
	NOTE='Compiling without M3D-C1 support'
	echo $NOTE | tee -a $LOG_FILE
else
	NOTE='Compiling with M3D-C1 support'
	echo $NOTE | tee -a $LOG_FILE
	echo "# $NOTE" >> include.mk
	echo "M3DC1_FLAG     = -DM3DC1" >> include.mk
	echo "M3DC1_INC      = -I $fusion_io_dir/include/_$fusion_io_arch" >> include.mk
        echo "M3DC1_LIBS     = -L $fusion_io_dir/lib/_$fusion_io_arch -lfusionio -lm3dc1 -lhdf5 -lstdc++" >> include.mk
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
echo "# Flags and libraries"				>> include.mk
echo 'FLAGS          = $(FGSL_CFLAGS) $(ODE_FLAGS)'	>> include.mk
echo 'LIBS           = $(FGSL_LIBS) $(M3DC1_LIBS)'	>> include.mk
echo ""							>> include.mk


#
echo 'FC             = $(COMPILER) $(FLAGS) -DFLARE $(OPT)' >> include.mk
echo 'FC_ADDONS      = $(COMPILER) $(FLAGS)         $(OPT)' >> include.mk
echo 'FC_DEBUG       = $(COMPILER) $(FLAGS) -DFLARE $(OPT_DEBUG)' >> include.mk


# generating header file config.h
echo "Base directory is" $HOME/$base_dir
echo "  character(*), parameter :: base_dir             =  '"$base_dir"'"    > config.h

echo "  character(*), parameter :: Boundary_input_file  =  'boundary.conf'" >> config.h
echo "  character(*), parameter :: Boundary_sub_dir     =  'boundary'"      >> config.h
echo "  character(*), parameter :: TAG                  =  'master'"        >> config.h
################################################################################
echo
