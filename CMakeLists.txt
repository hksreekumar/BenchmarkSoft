
############################################################################
# BenchmarkSoft CMAKE Module                                               #
############################################################################
#                                                                          #
# Author: Harikrishnan K. Sreekumar                                        #
#         Insitute for Acoustics                                           #
#         Technical University of Braunschweig                             #
#                                                                          #
# Description: Support for elPaSo Development                              #
############################################################################

project(BenchmarkSoft)

############################################################################
# 0. CMAKE Basics
############################################################################
cmake_minimum_required (VERSION 3.2)
cmake_policy(VERSION 3.2)
set( CMAKE_VERBOSE_MAKEFILE off )
############################################################################
# 1. Options
############################################################################
SET(LIB_HOME C:/software/libs)
############################################################################
# 1.1 Compiler Options
############################################################################
option(USE_GCC_COMPILER "Use GCC Compiler" ON )
set(GCC_VER 8)    # supports 4.8 , 7, 8
############################################################################
# 1.2 Math Libraries
############################################################################
option(USE_INTEL_MKL        "Use Intel MKL" ON)
option(USE_MPI              "Use MPI" OFF)
option(USE_PETSC            "Use PetSc" OFF)
set(PETSC_VER 3.11.4)    # Supports 3.7.7 , 3.8.4 , 3.9.4, 3.10.4
set(OMPI_VER 4.0.0)     # Supports 1.8.5 , 2.1.6 , 3.1.3, 4.0.0

add_definitions(-DPETSC_VER)
############################################################################
# 2. Setting C and C++ Compiler
############################################################################
IF("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
	IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64")
		SET(BMSOFT_OS_LIN_x86_64 FOUND)
	ENDIF()
ELSEIF("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
	IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "AMD64")
		SET(BMSOFT_OS_WIN_x86_64 FOUND)
	ENDIF()
ENDIF()

MESSAGE("> Build OS set to: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")

# ---- Set default compiler ------------------------------------------------#
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
IF(${BMSOFT_OS_LIN_x86_64})
	IF(${USE_GCC_COMPILER})
		set(CMAKE_C_COMPILER /usr/bin/gcc-${GCC_VER})
		set(CMAKE_CXX_COMPILER /usr/bin/g++-${GCC_VER})
	ENDIF()
ENDIF()

# ---- Identify target Compiler --------------------------------------------#
MESSAGE("> Identified compilers")
MESSAGE("              C   compiler | ID: ${CMAKE_C_COMPILER_ID} | VERSION: ${CMAKE_C_COMPILER_VERSION} | Path: ${CMAKE_C_COMPILER}")
MESSAGE("              CXX compiler | ID: ${CMAKE_CXX_COMPILER_ID} | VERSION: ${CMAKE_CXX_COMPILER_VERSION} | Path: ${CMAKE_CXX_COMPILER}")

set(BenchmarkSoft_SOURCES src/main.cpp)
MESSAGE("C compiler is ${CMAKE_C_COMPILER_ID} and has version ${CMAKE_C_COMPILER_VERSION} \nC++ compiler is ${CMAKE_CXX_COMPILER_ID} and has version ${CMAKE_CXX_COMPILER_VERSION} " ) 

add_definitions(-D_C_ID=\"${CMAKE_C_COMPILER_ID}\")
add_definitions(-D_C_VER=\"${CMAKE_C_COMPILER_VERSION}\")
add_definitions(-D_CXX_ID=\"${CMAKE_CXX_COMPILER_ID}\")
add_definitions(-D_CXX_VER=\"${CMAKE_CXX_COMPILER_VERSION}\")
add_definitions(-DUSE_GCC_COMPILER)
############################################################################
# 3. Setting Math Libraries
############################################################################
############################################################################
## 3.1 Setting Intel MKL
############################################################################
MESSAGE("Current source dir is " ${BenchmarkSoft_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${BenchmarkSoft_SOURCE_DIR}/cmake/modules)

IF (${USE_INTEL_MKL})
    find_package(MKL REQUIRED)
    add_definitions(-DUSE_INTEL_MKL)
    MESSAGE("MKL status: ${MKL_FOUND}")
    set(BenchmarkSoft_LIBS ${BenchmarkSoft_LIBS} ${MKL_LIBRARIES})    
ENDIF()
############################################################################
## 3.2 Setting OpenMPI
############################################################################
IF (${USE_MPI})
    find_package(MPI REQUIRED)
    add_definitions(-DUSE_MPI)
    MESSAGE("OMPI status: ${MPI_FOUND}")
    set(BenchmarkSoft_LIBS ${BenchmarkSoft_LIBS} ${MPI_LIBRARIES})
ENDIF()
############################################################################
## 3.3 Setting Petsc
############################################################################
IF (${USE_PETSC})
    find_package(PETSC REQUIRED)
    add_definitions(-DUSE_PETSC)
    MESSAGE("PETSC status: ${PETSC_FOUND}")
    set(BenchmarkSoft_LIBS ${BenchmarkSoft_LIBS} ${PETSC_LIBRARIES})
ENDIF()

############################################################################
# 4. SubDirectory CMAKE Linking
############################################################################
include(MacroAppendGlobalVariable)
add_subdirectory(src)
get_property(BenchmarkSoft_SOURCES GLOBAL PROPERTY BenchmarkSoft_SOURCES) 
get_property(BenchmarkSoft_INCLUDES GLOBAL PROPERTY BenchmarkSoft_INCLUDES) 

include_directories(${BenchmarkSoft_INCLUDES})
MESSAGE("BenchmarkSoft_LIBS: ${BenchmarkSoft_LIBS}")
MESSAGE("BenchmarkSoft_INCLUDES: ${BenchmarkSoft_INCLUDES}")

############################################################################
# 5. Executable Properties
############################################################################
add_executable(BenchmarkSoft ${BenchmarkSoft_SOURCES})
target_link_libraries(BenchmarkSoft ${BenchmarkSoft_LIBS})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin64)
