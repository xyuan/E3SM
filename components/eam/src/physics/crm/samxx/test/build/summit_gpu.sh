#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc/9.1.0 cuda/11.0.3 netcdf-c/4.8.0 netcdf-fortran/4.4.5 cmake/3.20.2  python/3.8-anaconda3 forge/20.1

unset ARCH
unset NCRMS
unset MACH
unset YAKL_DEBUG

export YAKL_DEBUG=true
export MACH="summit"
export NCHOME=${OLCF_NETCDF_C_ROOT}
export NFHOME=${OLCF_NETCDF_FORTRAN_ROOT}
export NCRMS=1
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS=" -g -O3 -ffree-line-length-none "
export CXXFLAGS=" -g -O3 "
export ARCH="CUDA"
export YAKL_CUDA_FLAGS="-arch sm_70 -g -O3  --use_fast_math -D__USE_CUDA__ -DKOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION -DTHRUST_IGNORE_CUB_VERSION_CHECK --expt-extended-lambda --expt-relaxed-constexpr"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_CUB_HOME="`pwd`/../../../../../../../../externals/cub"


