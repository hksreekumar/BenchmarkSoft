############################################################################
# BenchmarkSoft CMAKE Module for Intel MKL                                 #
############################################################################
#                                                                          #
# Author: Harikrishnan K. Sreekumar                                        #
#         Insitute for Acoustics                                           #
#         Technical University of Braunschweig                             #
#                                                                          #
# Description: Support for elPaSo Development                              #
############################################################################

MESSAGE("Finding MKL...")

set(MKLBASE_DIR "${LIB_HOME}")
set(MKLROOT_DIR "${MKLBASE_DIR}/intel/mkl")
set(MKLROOT "${MKLBASE_DIR}/intel/mkl")
message("MKLROOT_DIR is: ${MKLROOT_DIR}")

set(MKL_INCLUDE_DIR "${MKLBASE_DIR}/intel/mkl/include")
message("MKL_INCLUDE_DIR is: ${MKL_INCLUDE_DIR}")

set(MKL_LIB_DIR "${MKLBASE_DIR}/intel/mkl/lib/intel64")
message("MKL_LIB_DIR is: ${MKL_LIB_DIR}")

#set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
find_library(MKL_INTEL_ILP64_LIBRARY libmkl_intel_ilp64.a ${MKL_LIB_DIR})
find_library(MKL_INTEL_THREAD_LIBRARY libmkl_intel_thread.a ${MKL_LIB_DIR})
find_library(MKL_CORE_LIBRARY libmkl_core.a ${MKL_LIB_DIR})
find_library(MKL_OMP_LIBRARY libiomp5.a ${MKLBASE_DIR}/intel/compilers_and_libraries/linux/lib/intel64)

set (MKL_LIBRARIES "-Wl,--start-group ${MKL_INTEL_ILP64_LIBRARY} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} -Wl,--end-group ${MKL_OMP_LIBRARY} -lpthread -lm -ldl")

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64 -m64 -I${MKLROOT_DIR}/include")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL REQUIRED_VARS MKL_INTEL_ILP64_LIBRARY MKL_INTEL_THREAD_LIBRARY MKL_CORE_LIBRARY MKL_OMP_LIBRARY)

MESSAGE("Finding MKL...Complete.")
MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
