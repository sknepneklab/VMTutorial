
CMAKE_MINIMUM_REQUIRED(VERSION 3.14 FATAL_ERROR)

project (VMTOOLKIT C CXX)

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
  ${CMAKE_SOURCE_DIR}/VMToolkit/VM/src/version.hpp.in
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
set(VMTOOLKIT_LIBS ${VTK_LIBS})


# ##############################################
# place all executables in the build directory 
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)
add_subdirectory(extern/pybind11)

pybind11_add_module(VMTOOLKIT ${VMTOOLKIT_SRCS})


if (SKBUILD)
    execute_process(COMMAND python -m site --user-site OUTPUT_VARIABLE PYTHON_SITE OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "Local Python site ${PYTHON_SITE}")
    set(lib_path "${PYTHON_SITE}/VMToolkit/VM")
else()
    set(lib_path "${CMAKE_BINARY_DIR}")
endif()

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${lib_path}" is_system)
if ("${is_system}" STREQUAL "-1")
    target_link_libraries(VMTOOLKIT PRIVATE ${VMTOOLKIT_LIBS})
    install(TARGETS VMTOOLKIT DESTINATION .)
    set_target_properties(VMTOOLKIT PROPERTIES PREFIX ""  OUTPUT_NAME "vm" SUFFIX ".so" CXX_STANDARD 14  INSTALL_RPATH_USE_LINK_PATH TRUE INSTALL_RPATH "${lib_path}")
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    if (VTKLocal)
      add_dependencies(VMTOOLKIT VTKLocal)
      install(DIRECTORY "${lib_path}/extern/VTK" DESTINATION extern)
  endif()
      
endif()

#set_target_properties(VMTOOLKIT PROPERTIES PREFIX ""  OUTPUT_NAME "vm" SUFFIX ".so" BUILD_RPATH_USE_ORIGIN TRUE CXX_STANDARD 14)
#set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
#if (VTKLocal)
#  add_dependencies(VMTOOLKIT VTKLocal)
#  install(DIRECTORY ${CMAKE_BINARY_DIR}/extern/VTK DESTINATION extern)
#endif()


