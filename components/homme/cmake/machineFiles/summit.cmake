#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash

#module load cmake/3.6.1 cuda/9.1.76 gcc/6.4.0 netlib-lapack/3.6.1
#module load netcdf/4.6.1 netcdf-fortran/4.4.4
#module load hdf5/1.10.3

SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{OLCF_HDF5_ROOT} CACHE FILEPATH "")

#does not help
#set(NetCDF_C_LIBRARY     $ENV{OLCF_NETCDF_ROOT}/lib )
#set(NetCDF_C_INCLUDE_DIR $ENV{OLCF_NETCDF_ROOT}/include )

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(ENABLE_CUDA FALSE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(USE_TRILINOS OFF CACHE BOOL "")
SET(KOKKOS_PATH "$ENV{HOME}/kokkos/build-serial-cuda-nodebug" CACHE STRING "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "/ccs/home/onguba/kokkos/bin/nvcc_wrapper" CACHE STRING "")

set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

set (USE_NUM_PROCS 4 CACHE STRING "")

#set (OPT_FLAGS "-mcpu=power9 -mtune=power9" CACHE STRING "")
SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")