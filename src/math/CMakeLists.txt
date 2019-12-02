#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp *.c)
MACRO_APPEND_GLOBAL_VARIABLE(BenchmarkSoft_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(BenchmarkSoft_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
FILE(GLOB BenchmarkSoft_src_sub_subdirectories "*")
FOREACH(sf ${BenchmarkSoft_src_sub_subdirectories})
  IF(NOT sf MATCHES "\\.txt" AND IS_DIRECTORY ${sf} )
    add_subdirectory(${sf})
  ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#
