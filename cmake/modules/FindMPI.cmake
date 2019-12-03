############################################################################
# BenchmarkSoft CMAKE Module for MPI                               #
############################################################################
#                                                                          #
# Author: Harikrishnan K. Sreekumar                                        #
#         Insitute for Acoustics                                           #
#         Technical University of Braunschweig                             #
#                                                                          #
# Description: Support for elPaSo Development                              #
############################################################################

MESSAGE("Finding MPI version ${OMPI_VER}...")

set(MPIROOT_DIR "/opt/BS/openmpi-${OMPI_VER}")
set(MPIROOT "/opt/BS/openmpi-${OMPI_VER}")
message("MPIROOT_DIR is: ${MPIROOT_DIR}")

set(MPI_INCLUDE_DIR "${MPIROOT_DIR}/gcc-opt/include")
message("MPI_INCLUDE_DIR is: ${MPI_INCLUDE_DIR}")

set(MPI_LIB_DIR "${MPIROOT_DIR}/gcc-opt/lib")
message("MPI_LIB_DIR is: ${MPI_LIB_DIR}")

find_library(LIB_MPI libmpi.so ${MPI_LIB_DIR} NO_DEFAULT_PATH )
#find_library(LIB_MPI_CXX libmpi_cxx.so ${MPI_LIB_DIR} NO_DEFAULT_PATH)

SET(MPI_LIBRARIES "${LIB_MPI}") 
#SET(MPI_LIBRARIES "${LIB_MPI};${LIB_MPI_CXX}") 
#SET(MPI_LIBRARIES "${LIB_MPI_CXX}") 

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI REQUIRED_VARS LIB_MPI)

MESSAGE("Finding MPI...Complete.")
MESSAGE("MPI_LIBRARIES: ${MPI_LIBRARIES}")
