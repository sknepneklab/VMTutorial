#* ***************************************************************************
# *
# *  Copyright (C) 2017 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of AJM (Active Junction Model) program.
# *
# *  AJM is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *  (at your option) any later version.
# *
# *  AJM is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# * ****************************************************************************/

CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

project (VMTUTORIAL C CXX)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")

configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.hpp.in
  ${CMAKE_BINARY_DIR}/generated/version.hpp
)

include_directories(${CMAKE_BINARY_DIR}/generated)

# Set default CFlags
include (CMakeCFlagsSetup.txt)
# Configure some source files, include directories, and create variables listing all source files
include (CMakeSRCSetup.txt)
# Configure VTK libraries
include (CMakeVTKSetup.txt)
# Find Boost libraries
include (CMakeBoostSetup.txt)
# Find submodules (used to pull pybind11)
include (CMakeSubmodule.txt)

################################
## Define common libraries used by every target in MEMBRANE
#set(BOOST_LIBS  ${Boost_PROGRAM_OPTIONS_LIBRARY} )
set(VMTUTORIAL_LIBS ${VTK_LIBS})


# ##############################################
# place all executables in the build directory 
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


add_subdirectory(src)
add_subdirectory(extern/pybind11)

pybind11_add_module(VMTUTORIAL ${VMTUTORIAL_SRCS})


target_link_libraries(VMTUTORIAL PRIVATE ${VMTUTORIAL_LIBS})
install(TARGETS VMTUTORIAL DESTINATION ${CMAKE_INSTALL_PREFIX})

set_target_properties(VMTUTORIAL PROPERTIES PREFIX ""  OUTPUT_NAME "vmtutorial" SUFFIX ".so" CXX_STANDARD 14)
if (VTKLocal)
  add_dependencies(VMTUTORIAL VTKLocal)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/extern/VTK DESTINATION extern)
endif()


