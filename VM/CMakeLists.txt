
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
install(TARGETS VMTUTORIAL DESTINATION vm)

set_target_properties(VMTUTORIAL PROPERTIES PREFIX ""  OUTPUT_NAME "vm" SUFFIX ".so" BUILD_RPATH_USE_ORIGIN TRUE CXX_STANDARD 14)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
if (VTKLocal)
  add_dependencies(VMTUTORIAL VTKLocal)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/extern/VTK DESTINATION vm/extern)
endif()


