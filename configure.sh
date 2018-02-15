#!/bin/bash


################################################################################
# set default values

prefix=$HOME/.local
bindir=bin
libdir=lib
pythondir=$libdir/python
datadir=$HOME/Database/Magnetic_Configuration
EXTERNAL_DIR=external

# MPI wrapper for compiler
if type mpifort >/dev/null 2>/dev/null; then
   mpifc=mpifort
elif type mpif90 >/dev/null 2>/dev/null; then
   mpifc=mpif90
fi

with_fgsl=false
with_fio=false
hdf5="hdf5-openmpi"
python=python
f2py=f2py2.7
################################################################################


################################################################################
# scanning arguments
for opt in "$@"; do
    var=${opt%%=*}
    val=${opt##*=}
    case "$var" in
    -h|--help)
        echo -e "Configuration script for FLARE\n"
        echo -e "usage: ./configure.sh [options]\n"
        echo "Available options:"
        echo "  -h                      print short help message"
        echo "  --help                  print longer version of help message"
        echo
        echo "  --prefix=<path>         prefix for installation directories"
     if [ "$var" == "--help" ]; then
	echo "  --bindir=<path>         directory for installing FLARE executables"
        echo "  --datadir=<path>        data directory for magnetic configurations"
	echo "  --libdir=<path>         directory for the FLARE library"
	echo "  --pythondir=<path>      directory for python modules"
     fi
	echo
	echo "  --with-fgsl             compile with FGSL support"
     if [ "$var" == "--help" ]; then
	echo "  fgsl_cflags=<cflags>    compiler flags for compiling with FGSL support"
	echo "  fgsl_libs=<libs>        libraries for compiling with FGSL support"
        echo "  gsl_version=<version>   version of GSL library"
        echo "  fgsl_version=<version>  version of FGSL library"
	echo
     fi
	echo "  --with-fio              compile with Fusion-IO support"
     if [ "$var" == "--help" ]; then
        echo "  hdf5=<name>             name of MPI compatible HDF5 library"
	echo "  hdf5_cflags=<cflags>    compiler flags for compiling with HDF5 support"
	echo "  hdf5_libs=<libs>        libraries for compiling with HDF5 support"
        echo "  fio_dir=<path>          installation directory for Fusion-IO library"
        echo
     fi
     if [ "$var" == "--help" ]; then
        echo "  mpifc=<compiler>        MPI Fortran wrapper compiler"
        echo "  cflags=<cflags>         Extra flags to give to the Fortran compiler"
        echo "  python=<command>        Python 2.7 command"
        echo "  f2py=<command>          Fortran to Python 2.7 interface generator"
     fi

	exit 0
        ;;
    --prefix)
        prefix=$val
        ;;
    --bindir)
        bindir=$val
        ;;
    --datadir)
        datadir=$val
        ;;
    --libdir)
        libdir=$val
        ;;
    --pythondir)
        pythondir=$val
        ;;
    --with-fgsl)
        with_fgsl=true
        ;;
    fgsl_cflags)
        fgsl_cflags=$val
        ;;
    fgsl_libs)
        fgsl_libs=$val
        ;;
    gsl_version)
        gsl_version=$val
        ;;
    fgsl_version)
        fgsl_version=$val
        ;;
    --with-fio)
        with_fio=true
        ;;
    hdf5)
        hdf5=$val
        ;;
    hdf5_cflags)
        hdf5_cflags=$val
        ;;
    hdf5_libs)
        hdf5_libs=$val
        ;;
    fio_dir)
        fio_dir=$val
        ;;
    --without-gui)
        gui=
        ;;
    mpifc)
        mpifc=$val
        ;;
    cflags)
        cflags=$val
        ;;
    python)
        python=$val
        ;;
    f2py)
        f2py=$val
        ;;
    *)
        echo "error: invalid argument ${var}!"
        echo "see ./configure --help"
        exit 1
    esac
done
echo

bindir=$prefix/$bindir
libdir=$prefix/$libdir
pythondir=$prefix/$pythondir
################################################################################


################################################################################
# check configuration
function check_command() {
   message="checking for $(basename $1)... "
   if ! type "$1" > /dev/null 2>&1; then
      echo $message "error: cannot find $1!"
      exit 1
   fi
   echo $message $1 | tee -a $LOG_FILE
}

# pkg-config
if type pkg-config >/dev/null 2>/dev/null; then
   pkg_config=pkg-config
fi

# set variable from output of pkg-config
function autoset() {
   # check if variable is already set
   eval val=\$$1
   if [ -n "$val" ]; then
      return
   fi

   # set option for pkg-config
   opt=$2
   if [ -z "$opt" ]; then
      echo "error in autoset: no option given for pkg-config!"
      exit 1
   fi

   # set package name
   pkg=$3
   if [ -z "$pkg" ]; then
      echo "error in autoset: no package given for pkg-config!"
      exit 1
   fi
   if ! type $pkg_config >/dev/null 2>/dev/null; then
      echo "error: pkg-config is not available for information about ${pkg}!"
      echo "a value for variable $1 must be given as option for configure.sh!"
      exit 1
   fi
   if ! $pkg_config --exists $pkg; then
      echo "error: package $pkg does not exist!"
      echo "a value for variable $1 must be given as option for configure.sh!"
      exit 1
   fi

   # set variable from output of pkg-config
   val=`$pkg_config --$opt $pkg`
   #echo "RUNNING $pkg_config --$opt $pkg: $val"
   eval $1='$val'
}


# initialize configuration and log files ---------------------------------------
incmk=include.mk;  rm -f $incmk;  touch $incmk
cnf=config.h;      rm -f $cnf;    touch $cnf

LOG_FILE=configure.log
echo "command: $0 $@" > $LOG_FILE
echo "timestamp: $(date)" >> $LOG_FILE

# Main program
echo "# Main program"                                          >> $incmk
echo "PROGRAM        = FLARE"                                  >> $incmk
echo "FLARELIB       = libFLARE.so"                            >> $incmk
echo >> $incmk
#-------------------------------------------------------------------------------


# MPI Fortran compiler ---------------------------------------------------------
if [ -z "$mpifc" ]; then
   echo "error: MPI Fortran compiler undefined, explicit use of argument mpifc required!"
   exit 1
fi
check_command $mpifc
# Fortran compiler to be used in MPI wrapper
FCLONG=$($mpifc -show | cut -d\  -f1)
compiler=$(basename $FCLONG)
echo "$mpifc is using $compiler" | tee -a $LOG_FILE
# set compiler flags (if necessary)
if [ -z "$cflags" ]; then
case "$compiler" in
gfortran)
   cflags="-O2 -fconvert=big-endian -fPIC"
   cflags_debug="-g -fconvert=big-endian -fPIC -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace"
   ;;
ifort)
   cflags="-O2 -convert big_endian -fPIC"
   cflags_debug="-g -convert big_endian -fPIC -check all -debug all -fp-model strict -fp-speculation=off -fp-stack-check -fpconstant -fpe0 -traceback -fstack-security-check"
   ;;
*)
   echo "error: cannot find flags for Fortran compiler $compiler, explicit user of cflags required!"
   exit 1
esac
fi

echo "# Fortran compiler and options"                          >> $incmk
echo "MPIFC          = $mpifc -DMPI $cflags"                   >> $incmk
echo "MPIFC_DEBUG    = $mpifc -DMPI $cflags_debug"             >> $incmk
echo >> $incmk
#-------------------------------------------------------------------------------


# GNU Scientific Library -------------------------------------------------------
if $with_fgsl; then
   NOTE="Compiling with FGSL support"
   echo $NOTE | tee -a $LOG_FILE

   autoset fgsl_cflags   cflags      fgsl
   autoset fgsl_libs     libs        fgsl
   autoset gsl_version   modversion  gsl
   autoset fgsl_version  modversion  fgsl

   gsl_version_major=$(echo $gsl_version | cut -d. -f1)
   gsl_version_minor=$(echo $gsl_version | cut -d. -f2)
   gsl_version_micro=$(echo $gsl_version | cut -d. -f3)

   fgsl_version_major=$(echo $fgsl_version | cut -d. -f1)
   fgsl_version_minor=$(echo $fgsl_version | cut -d. -f2)
   fgsl_version_micro=$(echo $fgsl_version | cut -d. -f3)

   echo "#define FGSL" >> config.h
   echo "#define GSL_VERSION_MAJOR_FORTRAN $gsl_version_major" >> $cnf
   echo "#define GSL_VERSION_MINOR_FORTRAN $gsl_version_minor" >> $cnf

   echo "# GNU Scientific Library"                             >> $incmk
   echo "FGSL_CFLAGS    = $fgsl_cflags"                        >> $incmk
   echo "FGSL_LIBS      = $fgsl_libs"                          >> $incmk
   echo >> $incmk
else
   NOTE="Compiling without FGSL support"
   echo $NOTE | tee -a $LOG_FILE
fi
#-------------------------------------------------------------------------------


# ODE solvers ------------------------------------------------------------------
echo "# additional ODE solvers" >> $incmk
declare -a ode_solver_files=("dlsode.f")
declare -a ode_solver_names=("DLSODE"  )
n_ode=${#ode_solver_files[@]}

for (( i=0; i<${n_ode}; i++ )); do
	if [ -f src/$EXTERNAL_DIR/${ode_solver_files[$i]} ]; then
		NOTE="Compiling with ${ode_solver_names[$i]}"
		echo $NOTE | tee -a $LOG_FILE
		ode_cflags="$ode_cflags -D${ode_solver_names[$i]}"
		ode_objects="$ode_objects ${ode_solver_files[$i]%%.f*}.o"
	else
		NOTE="Compiling without ${ode_solver_names[$i]}"
		echo $NOTE | tee -a $LOG_FILE
	fi
done
echo "ODE_CFLAGS     = $ode_cflags"                            >> $incmk
echo "ODE_OBJECTS    = $ode_objects"                           >> $incmk
echo >> $incmk
#-------------------------------------------------------------------------------


# Fusion-IO installation -------------------------------------------------------
if $with_fio; then
   NOTE='Compiling with Fusion-IO support'
   echo $NOTE | tee -a $LOG_FILE
   echo "# $NOTE"                                              >> $incmk

   # check HDF5 installation
   autoset hdf5_cflags  cflags  $hdf5
   autoset hdf5_libs    libs    $hdf5
   echo "HDF5_CFLAGS    = $hdf5_cflags"                        >> $incmk
   echo "HDF5_LIBS      = $hdf5_libs"                          >> $incmk

   # check for fusionio library
   if [ -z "$fio_dir" ]; then
      if ! ldconfig -p | grep fusionio; then
         echo "error: fusionio library not found!"
         exit 1
      fi
   else
      fio_so=$fio_dir/lib/libfusionio.so
      fio_a=$fio_dir/lib/libfusionio.a
      if [ ! -f $fio_so ] && [ ! -f $fio_a ]; then
         echo "error: cannot find fusionio library in ${fio_dir}/lib!"
         exit 1
      fi
   fi
   fio_cflags="-DFIO"
   fio_libs="-lfusionio -lm3dc1 -lstdc++"
   if [ -n "$fio_dir" ]; then
      fio_cflags="$fio_cflags -I${fio_dir}/include"
      fio_libs="$fio_libs -L${fio_dir}/lib"
      fio_path="-Wl,-rpath-link,${fio_dir}/lib"
   fi

   echo "FIO_CFLAGS     = $fio_cflags"                         >> $incmk
   echo "FIO_LIBS       = $fio_libs"                           >> $incmk
   echo "FIO_PATH       = $fio_path"                           >> $incmk
   echo >> $incmk
else
   NOTE='Compiling without Fusion-IO support'
   echo $NOTE | tee -a $LOG_FILE
fi
#-------------------------------------------------------------------------------


# f2py--------------------------------------------------------------------------
check_command $python
check_command $f2py
echo "# f2py"                                                  >> $incmk
echo "PYTHON         = $python"                                >> $incmk
echo "F2PY           = $f2py --fcompiler=$compiler"            >> $incmk
echo >> $incmk
#-------------------------------------------------------------------------------
################################################################################


################################################################################
# installation directories
echo "# set up directories for executables and data"           >> $incmk
echo "BINDIR         = $bindir"                                >> $incmk
echo "LIBDIR         = $libdir"                                >> $incmk
echo "PYTHONDIR      = $pythondir"                             >> $incmk
echo "Binary directory is $bindir" | tee -a $LOG_FILE
echo "Library directory is $libdir" | tee -a $LOG_FILE
echo "Python module directory is $pythondir" | tee -a $LOG_FILE

# data directory
echo "#define DATABASE_DIR '$datadir'" >> $cnf
echo "DATADIR        = $datadir"                               >> $incmk
echo "Data directory is $datadir" | tee -a $LOG_FILE
echo >> $incmk
echo
#-------------------------------------------------------------------------------


# set up local source directories ----------------------------------------------
echo "# set up local source directories"                       >> $incmk
echo "EXTERNAL_DIR   = $EXTERNAL_DIR"                          >> $incmk
echo "CORE_DIR       = core"                                   >> $incmk
echo "BFIELD_DIR     = bfield"                                 >> $incmk
echo "GEOMETRY_DIR   = geometry"                               >> $incmk
echo "GRIDGEN_DIR    = fieldline_grid"                         >> $incmk
echo "TOOLS_DIR      = tools"                                  >> $incmk
echo "DEVEL_DIR      = development"                            >> $incmk
echo "ADDONS_DIR     = addons"                                 >> $incmk
echo >> $incmk
#-------------------------------------------------------------------------------


# Flags and libraries ----------------------------------------------------------
CFLAGS='$(FGSL_CFLAGS) $(HDF5_CFLAGS) $(FIO_CFLAGS) $(ODE_CFLAGS)'
LIBS='$(FGSL_LIBS)   $(HDF5_LIBS)   $(FIO_LIBS)'

echo "# Flags and libraries"                                   >> $incmk
echo "CFLAGS         = $CFLAGS"                                >> $incmk
echo "LIBS           = $LIBS"                                  >> $incmk
echo 'FC             = $(MPIFC)       $(CFLAGS) -DFLARE'       >> $incmk
echo 'FC_ADDONS      = $(MPIFC)       $(CFLAGS)'               >> $incmk
echo 'FC_DEBUG       = $(MPIFC_DEBUG) $(CFLAGS) -DFLARE'       >> $incmk
################################################################################
