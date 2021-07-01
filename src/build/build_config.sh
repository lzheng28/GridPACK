rm -rf CMake*
CC=gcc
CXX=g++
CXXFLAGS=-lrt
export CC CXX CXXFLAGS

cmake \
    -D PETSC_ARCH:STRING="arch-linux2-cxx-debug" \
    -D PETSC_DIR:STRING="/home/lei/Install_package/petsc-v3.7.7" \
    -D GA_DIR:STRING='/home/lei/software/ga-5.7' \
    -D PARMETIS_DIR:PATH="/usr" \
    -D MPI_CXX_COMPILER:STRING="mpicxx" \
    -D MPI_C_COMPILER:STRING="mpicc" \
    -D MPIEXEC:STRING="mpiexec" \
    -D MPIEXEC_MAX_NUMPROCS:STRING="2" \
    -D GRIDPACK_TEST_TIMEOUT:STRING=60 \
    -D USE_GLPK:BOOL=ON \
    -D GLPK_ROOT_DIR:PATH="/usr" \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D CMAKE_INSTALL_PREFIX:PATH="$HOME/software/gridpack_lzheng" \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D BOOST_ROOT:STRING='/home/lei/software/boost_1_70_0' \
    -D HELICS_INSTALL_DIR:PATH="/home/lei/software/helics-2.6.1" \
    -D ZEROMQ_INSTALL_DIR:PATH="/usr/lib/x86_64-linux-gnu" \
    CXXFLAGS=-lrt \
    ..
