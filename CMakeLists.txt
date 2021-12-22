cmake_minimum_required(VERSION 3.0.1)

project(AsteroidProject
        VERSION 0.0)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 CXX_STANDARD_14_SUPPORT -g)
if(CXX_STANDARD_14_SUPPORT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
endif()


option(CORAX_TESTS_ENABLED "Enable and build the tests" OFF)
set(CORAX_TESTS_ENABLED OFF)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
add_subdirectory(coraxlib)

if (NOT DISABLE_MPI)
  include(FindMPI)
  if(MPI_FOUND)
    set(WITH_MPI TRUE)
    add_definitions(-DWITH_MPI)
    message(STATUS "WITH MPI ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()


set(ASTEROID_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/core
   ${CMAKE_CURRENT_SOURCE_DIR}/coraxlib
   ${MPI_INCLUDE_PATH}
   ${CORAX_INCLUDE_DIRS})
  

add_subdirectory(src)

