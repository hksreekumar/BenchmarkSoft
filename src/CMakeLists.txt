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

#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(BenchmarkSoft_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(BenchmarkSoft_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
FILE(GLOB BenchmarkSoft_src_subdirectories "*")
FOREACH(f ${BenchmarkSoft_src_subdirectories})
  IF(NOT f MATCHES "\\.txt" AND IS_DIRECTORY ${f} )
    add_subdirectory(${f})
  ENDIF()
ENDFOREACH()

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

get_property(BenchmarkSoft_INCLUDES GLOBAL PROPERTY BenchmarkSoft_INCLUDES) 
get_property(BenchmarkSoft_SOURCES GLOBAL PROPERTY BenchmarkSoft_SOURCES) 

include_directories(${BenchmarkSoft_INCLUDES})
LIST(REMOVE_ITEM BenchmarkSoft_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

MESSAGE("BenchmarkSoft_INCLUDES: ${BenchmarkSoft_INCLUDES}")
MESSAGE("BenchmarkSoft_SOURCES: ${BenchmarkSoft_SOURCES}")

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()
