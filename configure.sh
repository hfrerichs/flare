#!/bin/bash


################################################################################
# predefine parameter
base_dir="Database/Magnetic_Configuration"

EMC3_dir=""

fusion_io_dir=$M3DC1_INSTALL_DIR
fusion_io_arch=`uname`
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
	echo "  --EMC3_dir=DIR          set EMC3 source directory"
	echo ""
	echo "  --fusion_io_dir=DIR     set directory of M3D-C1 / fusion_io installation"
	echo "  --fusion_io_arch=ARCH   set architecture string used for M3D-C1 (if not equal to `uname`)"
	echo ""

	exit
    elif [ "$par" == "--base_dir" ]; then
        base_dir=$val
    elif [ "$par" == "--EMC3_dir" ]; then
        EMC3_dir=$val
    elif [ "$par" == "--fusion_io_dir" ]; then
        fusion_io_dir=$val
    elif [ "$par" == "--fusion_io_arch" ]; then
        fusion_io_arch=$val
    else
        echo "error: unkown parameter " $par
    fi
done
################################################################################


################################################################################
# generate configuration file
rm -rf include.mk
touch include.mk

LOG_FILE=configure.log
echo "timestamp: $(date)" > $LOG_FILE

# MPI support / Fortran compiler
if type mpif90 >/dev/null 2>/dev/null; then
	echo "Compiling with MPI support" | tee -a $LOG_FILE
	echo "FC_BASE        = mpif90 -DparallelMPI" >> include.mk
elif type ifort >/dev/null 2>/dev/null; then
	echo "Using Inter Fortran compiler" | tee -a $LOG_FILE
	echo "FC_BASE        = ifort" >> include.mk
elif type gfortran >/dev/null 2>/dev/null; then
	echo "Using GNU Fortran compiler" | tee -a $LOG_FILE
	echo "FC_BASE        = gfortran" >> include.mk
fi
echo "OPT            = -O2" >> include.mk
echo "OPT_DEBUG      = -g" >> include.mk
echo "" >> include.mk


# setting local source directories
echo "# setting local source directories" >> include.mk
echo "EXTERNAL_DIR   = external" >> include.mk
echo "EMC3_LINK_DIR  = \$(EXTERNAL_DIR)/EMC3" >> include.mk
echo "MATH_DIR       = math" >> include.mk
echo "GEOMETRY_DIR   = geometry" >> include.mk
echo "BFIELD_DIR     = bfield" >> include.mk
echo "PFC_DIR        = pfc" >> include.mk
echo "TOOLS_DIR      = tools" >> include.mk
echo "" >> include.mk


# checking for coupling to EMC3
if [ "$EMC3_dir" == "" ]; then
	NOTE='Compiling without support for field aligned grid input (bases on EMC3 sources)'
	echo $NOTE | tee -a $LOG_FILE
else
	NOTE='Compiling with support for field aligned grid input (bases on EMC3 sources)'
	echo $NOTE | tee -a $LOG_FILE
	echo "# $NOTE" >> include.mk
	echo "EMC3_FLAG      = -DEMC3" >> include.mk
        echo "EMC3_SRC_DIR   = $EMC3_dir/MAIN" >> include.mk
	if [ ! -d "$EMC3_dir/MAIN" ]; then
		echo "error: EMC3 source files not found!"
		exit
	fi
	echo "EMC3_OBJ       = \
               PHYS_CONST.o\
               GEOMETRY_PL.o\
               SURFACE_PL.o\
               MAPPING.o\
               CELL_GEO.o\
               cell_def.o\
               cell_user.o\
               check.o\
               ibm_iface.o\
               random.o\
               real_to_ft.o\
               service.o\
               sf_def_user.o\
               sf_jump.o" >> include.mk
	echo "EMC3_OBJ_LONG  = \$(addprefix \$(EMC3_LINK_DIR)/,\$(EMC3_OBJ))" >> include.mk
	echo "" >> include.mk


fi


# checking for M3D-C1/fusion_io installation
if [ "$fusion_io_dir" == "" ]; then
	NOTE='Compiling without M3D-C1 support'
	echo $NOTE | tee -a $LOG_FILE
else
	NOTE='Compiling with M3D-C1 support'
	echo $NOTE | tee -a $LOG_FILE
	echo $NOTE
	echo "# $NOTE" >> include.mk
	echo "M3DC1_FLAG     = -D_M3DC1" >> include.mk
	echo "M3DC1_INC      = -I $fusion_io_dir/include/_$fusion_io_arch" >> include.mk
        echo "LIBS           = -L $fusion_io_dir/lib/_$fusion_io_arch -lfusionio -lm3dc1 -lhdf5 -lstdc++" >> include.mk
	echo "" >> include.mk
fi


echo 'FC             = $(FC_BASE) -DFLARE $(OPT)' >> include.mk
echo 'FC_DEBUG       = $(FC_BASE) -DFLARE $(OPT_DEBUG)' >> include.mk
echo "Setting base directory to" '$HOME'/$base_dir
echo "      character(*), parameter :: base_dir = '"$base_dir"'" > config.h
################################################################################
