cmake_minimum_required(VERSION 3.14)  # Specify your required version
project(VMTutorial LANGUAGES CXX)

# Optional: Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Add VM subdirectory
add_subdirectory(VMToolkit/VM)