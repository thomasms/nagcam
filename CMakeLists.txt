# CMake project file for forts
cmake_minimum_required (VERSION 2.8)

# project name and details
project (NagCam)
enable_language (Fortran)
set(PROJECT_DESCRIPTION "NAGCAM")
set(PROJECT_URL "https://git.ccfe.ac.uk/thomasms/nagcam")

# set source directories
set (PROJ_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
set (PROJ_BUILD_DIR ${CMAKE_BINARY_DIR})

# library names
set (PROJ_LIB_NAME nagcamlib)

# Set release build by default
set (CMAKE_BUILD_TYPE Release CACHE STRING "Select which configuration to build.")

# set the general fortran Debug/Release compiler flags
set(PROJ_RELEASE_F90_FLAGS "-O3")
set(PROJ_RELEASE_CXX_FLAGS "-O3")
set(PROJ_DEBUG_F90_FLAGS "-g -O0")
set(PROJ_DEBUG_CXX_FLAGS "-g -O0")

set(CMAKE_Fortran_FLAGS_RELEASE ${PROJ_RELEASE_FLAGS})
set(CMAKE_Fortran_FLAGS_DEBUG ${PROJ_DEBUG_FLAGS})

# C++ flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fvisibility=hidden -fvisibility-inlines-hidden")
set(CMAKE_CXX_FLAGS_DEBUG ${PROJ_DEBUG_CXX_FLAGS})
set(CMAKE_CXX_FLAGS_RELEASE ${PROJ_RELEASE_CXX_FLAGS})

# Fortran OS specific compiler flags
if (CMAKE_Fortran_COMPILER_ID MATCHES MSVC)
    set ( PROJ_DEBUG_F90_FLAGS "${PROJ_DEBUG_F90_FLAGS} -fp-model")
    set ( PROJ_RELEASE_FLAGS  "${PROJ_RELEASE_F90_FLAGS} -fp-model" )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set ( PROJ_DEBUG_F90_FLAGS "${PROJ_DEBUG_F90_FLAGS} -debug  -CB -CA -CU -std08 -fpp -fp-model source")
    set ( PROJ_RELEASE_FLAGS  "${PROJ_RELEASE_F90_FLAGS} -std08 -unroll -fpp -fp-model source -diag-disable 8291,8577" )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set ( PROJ_DEBUG_F90_FLAGS "${PROJ_DEBUG_F90_FLAGS} -cpp -ffree-form  -ffloat-store -fno-sign-zero -std=f2008 -ffree-line-length-none")
    set ( PROJ_RELEASE_FLAGS  "${PROJ_RELEASE_F90_FLAGS} -cpp -ffree-form -ffloat-store -fno-sign-zero -std=f2008 -ffree-line-length-none" )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES NAG)
    set ( PROJ_DEBUG_F90_FLAGS "${PROJ_DEBUG_F90_FLAGS} -fpp -free -float-store -f2008")
    set ( PROJ_RELEASE_FLAGS  "${PROJ_RELEASE_F90_FLAGS} -fpp -free -float-store -f2008")
else ()
    message(ERROR "OS >>> Fortran Compiler not supported!!")
endif ()

# Release or Debug build
if(CMAKE_BUILD_TYPE MATCHES Release)
    message(STATUS "Release build")
    set(F90_FLAGS ${PROJ_RELEASE_F90_FLAGS})
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Debug build")
    set(F90_FLAGS ${PROJ_DEBUG_F90_FLAGS})
else()
    message(STATUS "Default build (Release)")
    set(F90_FLAGS ${PROJ_RELEASE_F90_FLAGS})
endif()

# set the output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# Collect all the mod files into their own include directory
set(PROJ_MODULE_DIR ${PROJ_BUILD_DIR}/include)

# set include and link directories
include_directories (${PROJ_MODULE_DIR})
link_directories (${PROJ_BUILD_DIR}/lib)

# sources
add_subdirectory(${PROJ_SOURCE_DIR})


