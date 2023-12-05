# * ***************************************************************************
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


macro(add_source_dir dir)
add_subdirectory(${dir})

# find all source files in this directory
file(GLOB VMTOOLKIT_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp)
file(GLOB VMTOOLKIT_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.hpp)

# add them to the source group
source_group(${dir} FILES ${VMTOOLKIT_SRCS} ${VMTOOLKIT_HDRS})

# add them to the main samos list
list(APPEND _vmtutorial_sources ${VMTOOLKIT_SRCS} ${VMTOOLKIT_HDRS})
endmacro(add_source_dir)
#################



file(GLOB _VMTOOLKIT_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file(GLOB _VMTOOLKIT_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

SET(VMTOOLKIT_HDRS ${_VMTOOLKIT_HDRS} PARENT_SCOPE)
SET(VMTOOLKIT_SRCS ${_VMTOOLKIT_SRCS} PARENT_SCOPE)




