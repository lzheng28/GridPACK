#! /bin/sh


# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [-s] [name]"

set -- `getopt rds $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi

build="RelWithDebInfo"
shared="FALSE"
for o in $*; do
    case $o in
        -d)
            build="Debug"
            shift
            ;;
        -r)
            build="Release"
            shift
            ;;
        -s)
            shared="ON"
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "$0: error: $o: unknown option" >&2
            echo $usage >&2
            exit 2
    esac
done

if [ $# -gt 0 ]; then
    host="$1"
else
    host=`uname -n`
fi

rm -rf CMakeCache.txt CMakeFiles

options="-Wdev --debug-trycompile"

# useful build types: Debug, Release, RelWithDebInfo
common_flags="\
        -D BUILD_SHARED_LIBS:BOOL=$shared \
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
"

if [ $host == "flophouse" ]; then

    # RHEL 7 with GNU 4.8 compilers
    # The following are installed as packages:
    #   boost.x86_64                         1.53.0-26.el7
    #   boost-openmpi.x86_64                 1.53.0-26.el7
    #   (lots of other boost packages)
    #   

    prefix="/net/flophouse/files0/perksoft/linux64"
    PATH="${prefix}/bin:${PATH}"
    export PATH

    CC="/usr/bin/gcc"
    export CC
    CXX="/usr/bin/g++"
    export CXX

    if [ "$shared"x = "ON"x ]; then
        pdir="/net/flophouse/files0/perksoft/petsc-3.8.4"
        parch="rhel7-gnu48-real-opt-shared"
    else
        # pdir="/net/flophouse/files0/perksoft/petsc-3.4.5"
        # parch="rhel7-real-c++-static"
        # parch="rhel7-real-c-static"
        pdir="/net/flophouse/files0/perksoft/petsc-3.8.4"
        parch="rhel7-gnu48-real-opt"
        parch="rhel7-gnu48-complex-opt"
        parch="rhel7-gnu48-complex-opt-c"
    fi

    cplexroot="/opt/ibm/ILOG/CPLEX_Studio1261"

    cmake3 -Wdev --debug-trycompile \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D GA_DIR:PATH="$prefix/ga-c++" \
        -D BOOST_ROOT:STRING="/usr" \
        -D PETSC_DIR:STRING="$pdir" \
        -D PETSC_ARCH:STRING="$parch" \
        -D MPI_CXX_COMPILER:STRING="mpicxx" \
        -D MPI_C_COMPILER:STRING="mpicc" \
        -D MPIEXEC:STRING="mpiexec" \
        -D USE_CPLEX:BOOL=OFF \
        -D CPLEX_ROOT_DIR:PATH="$cplexroot" \
        -D USE_GLPK:BOOL=ON \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..

elif [ $host == "briareus" ]; then

    # Using GNU 4.9 and OpenMPI modules

    prefix="/files0/gridpack.gnu"
    PATH="${prefix}/bin:${PATH}"
    export PATH

    CC="/share/apps/gcc/4.9.2/bin/gcc"
    export CC
    CXX="/share/apps/gcc/4.9.2/bin/g++"
    export CXX

    cmake -Wdev --debug-trycompile \
        -D USE_PROGRESS_RANKS:BOOL=OFF \
        -D GA_DIR:PATH="$prefix" \
        -D BOOST_ROOT:STRING="$prefix" \
        -D PETSC_DIR:STRING="/files0/petsc-3.7.5" \
        -D PETSC_ARCH:STRING="gridpack-gnu-openmpi-real" \
        -D MPI_CXX_COMPILER:STRING="mpicxx" \
        -D MPI_C_COMPILER:STRING="mpicc" \
        -D MPIEXEC:STRING="mpiexec" \
        -D USE_CPLEX:BOOL=OFF \
        -D USE_GLPK:BOOL=OFF \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix" \
        $common_flags ..

elif [ $host == "WE32673" ]; then

    # Mac using CLang 6.0 compilers and MPICH via MacPorts
    # The following MacPorts packages are installed:
    #   clang-6.0 @6.0.1_0+analyzer+libstdcxx
    #   mpich-clang60 @3.2.1_4+gcc7
    #   boost @1.66.0_3+clang60+mpich+no_single+no_static+python27
    #   glpk @4.65_0
    #   doxygen @1.8.13_2+qt4+wizard
    # Global Arrays 5.7 built by hand
    # PETSc 3.8.4 w/ ParMETIS, SuperLU, etc., built by hand
    # Need to make sure the compiler set and MPI are selected, i.e.
    #   sudo port select clang mp-clang-6.0
    #   sudo port select mpi mpich-clang60-fortran
    # Cannot use PETSc < 3.8.0

    CC=/opt/local/bin/clang
    export CC
    CXX=/opt/local/bin/clang++
    export CXX

    prefix="/Users/d3g096/Projects/GridPACK"

    cmake $options \
        -D GA_DIR:STRING="$prefix" \
        -D BOOST_ROOT:STRING="/opt/local" \
        -D PETSC_DIR:PATH="$prefix/petsc-3.8.4" \
        -D PETSC_ARCH:STRING="arch-macosx-clang-real-shared-c" \
        -D MPI_CXX_COMPILER:STRING='/opt/local/bin/mpicxx' \
        -D MPI_C_COMPILER:STRING='/opt/local/bin/mpicc' \
        -D MPIEXEC:STRING='/opt/local/bin/mpiexec' \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
        -D USE_CPLEX:BOOL=OFF \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/opt/local" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/gridpack" \
        $common_flags ..

elif [ $host == "olympus.local" ]; then

    prefix="/pic/projects/gridpack/software"
    cmake $options \
        -D GA_DIR:STRING="/pic/projects/gridpack/ga-5-2" \
        -D GA_EXTRA_LIBS:STRING="-libverbs" \
	-D BOOST_ROOT:STRING="$prefix" \
	-D PETSC_DIR:STRING="$prefix/petsc-3.4.0" \
	-D PETSC_ARCH:STRING='olympus-openmpi-gnu-cxx-complex-opt' \
	-D MPI_CXX_COMPILER:STRING='mpicxx' \
	-D MPI_C_COMPILER:STRING='mpicc' \
	-D MPIEXEC:STRING='mpiexec' \
	$common_flags ..

elif [ $host == "tlaloc" ]; then

    # RHEL 6 with GNU 4.4 compilers w/ OpenMPI (available via EPEL)
    # Boost 1.55, GA 5.6.2, and PETSc 3.6.4 built by hand.  Available
    # Boost, GA, and PETSc packages were either too old or installed
    # in such a way as to be unrecognizable.

    prefix="/file0/perksoft"
    if [ "$shared"x = "ON"x ]; then
        pdir="/net/flophouse/files0/perksoft/petsc-3.8.4"
        parch="rhel6-complex-c-shared"
    else
        pdir="/net/flophouse/files0/perksoft/petsc-3.4.5"
        # pdir="/net/flophouse/files0/perksoft/petsc-3.4.2"
        pdir="/net/flophouse/files0/perksoft/petsc-3.5.4"
        pdir="/net/flophouse/files0/perksoft/petsc-3.8.4"
        parch="rhel6-real-c-static"
        #parch="rhel6-real-c++-static"
        parch="rhel6-complex-c-static"
        # parch="rhel6-complex-c++-static"
    fi
    
    cmake3 $options \
          -D GA_DIR:PATH="${prefix}/ga-c++" \
          -D BOOST_ROOT:PATH="${prefix}" \
          -D USE_PROGRESS_RANKS:BOOL=OFF \
          -D PETSC_DIR:PATH="$pdir" \
          -D PETSC_ARCH:STRING="$parch" \
          -D MPI_CXX_COMPILER:STRING="mpicxx" \
          -D MPI_C_COMPILER:STRING="mpicc" \
          -D MPIEXEC:STRING="mpiexec" \
          -D USE_GLPK:BOOL=OFF \
          -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
          -D GRIDPACK_TEST_TIMEOUT:STRING=10 \
          -D CMAKE_INSTALL_PREFIX:PATH="${prefix}/gridpack" \
          $common_flags ..

    

elif [ $host == "gridpackvm" ]; then

    prefix="$HOME/gridpack"

    CC=gcc
    CXX=g++
    export CC CXX

    cmake $options \
	-D PETSC_DIR:STRING="/usr/lib/petsc" \
	-D PARMETIS_DIR:PATH="/usr" \
	-D GA_EXTRA_LIBS:STRING="-lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -llapack -lblas -lgfortran" \
	-D MPI_CXX_COMPILER:STRING="mpicxx" \
	-D MPI_C_COMPILER:STRING="mpicc" \
	-D MPIEXEC:STRING="mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=20 \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/usr" \
        -D CMAKE_INSTALL_PREFIX:PATH="/usr" \
	-D HELICS_INSTALL_DIR:PATH="/home/huan495/HELICS/install-test" \
	$common_flags ..

elif [ $host == "debianvm" ]; then

    prefix="$HOME/gridpack"

    CC=gcc
    CXX=g++
    CFLAGS=-pthread
    CXXFLAGS=-pthread
    export CC CXX CFLAGS CXXFLAGS

    cmake $options \
	-D PETSC_DIR:STRING="/usr/lib/petsc" \
	-D PARMETIS_DIR:PATH="/usr" \
	-D GA_EXTRA_LIBS:STRING="-lscalapack-openmpi -lblacs-openmpi -llapack -lblas -lgfortran" \
	-D MPI_CXX_COMPILER:STRING="mpicxx" \
	-D MPI_C_COMPILER:STRING="mpicc" \
	-D MPIEXEC:STRING="mpiexec" \
        -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
        -D GRIDPACK_TEST_TIMEOUT:STRING=30 \
        -D USE_GLPK:BOOL=ON \
        -D GLPK_ROOT_DIR:PATH="/usr" \
        -D BUILD_SHARED_LIBS:BOOL=OFF \
        -D CMAKE_INSTALL_PREFIX:PATH="/usr" \
	$common_flags ..

else

    echo "Unknown host: $host"
    exit 2
    
fi

