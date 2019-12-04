############################################################################
# BenchmarkSoft CMAKE Module for PETSC                                     #
############################################################################
#                                                                          #
# Author: Harikrishnan K. Sreekumar                                        #
#         Insitute for Acoustics                                           #
#         Technical University of Braunschweig                             #
#                                                                          #
# Description: Support for elPaSo Development                              #
############################################################################

MESSAGE("Finding PETSC version ${PETSC_VER}...")

set(PETSCROOT_DIR "${LIB_HOME}/petsc-${PETSC_VER}")
set(PETSCROOT "${LIB_HOME}/petsc-${PETSC_VER}")
message("PETSCROOT_DIR is: ${PETSCROOT_DIR}")

set(PETSC_INCLUDE_DIR "${PETSCROOT_DIR}/include;${PETSCROOT_DIR}/gcc-cxx-complex-o/include")
message("PETSC_INCLUDE_DIR is: ${PETSC_INCLUDE_DIR}")

set(PETSC_LIB_DIR "${PETSCROOT_DIR}/gcc-cxx-complex-o/lib")
message("PETSC_LIB_DIR is: ${PETSC_LIB_DIR}")

find_library(LIB_PETSC libpetsc.so ${PETSC_LIB_DIR})

SET(PETSC_LIBRARIES ${LIB_PETSC}) 

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PETSCROOT_DIR}/include -I${PETSCROOT_DIR}/gcc-cxx-complex-o/include")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC REQUIRED_VARS LIB_PETSC)

MESSAGE("Finding PETSC...Complete.")
MESSAGE("PETSC_LIBRARIES: ${PETSC_LIBRARIES}")
