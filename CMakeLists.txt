# 
#  Copyright (C) 2017 cogimon
#  Author: Enrico Mingo Hoffman
#  email: enrico.mingo@iit.it
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program. If not, see <http://www.gnu.org/licenses/>
#

cmake_minimum_required(VERSION 2.8.12)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 OLD)
    cmake_policy(SET CMP0005 NEW)
    cmake_policy(SET CMP0017 NEW)
endif(COMMAND cmake_policy)
include(ExternalProject)
PROJECT(compliant_stabilizer)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(MacroInstallLib)

FIND_PACKAGE(YCM REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(OpenSoT QUIET)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wmaybe-uninitialized -Wuninitialized")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wmaybe-uninitialized -Wuninitialized")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

INCLUDE_DIRECTORIES(include ${EIGEN3_INCLUDE_DIR})

set(compliant_stabilizer_sources src/compliantstabilizer.cpp
                          src/FilterHClass.cpp)

if (${OpenSoT_FOUND})
    message("Adding src/CoMStabilizer.cpp to compilation")
    set(compliant_stabilizer_sources ${compliant_stabilizer_sources}
        src/CoMStabilizer.cpp)

    INCLUDE_DIRECTORIES(${OpenSoT_INCLUDE_DIRS})

endif()

ADD_LIBRARY(compliant_stabilizer SHARED ${compliant_stabilizer_sources})

if (${OpenSoT_FOUND})
    TARGET_LINK_LIBRARIES(compliant_stabilizer ${OpenSoT_LIBRARIES})
endif()

## Addons for libCartesianOpenSot
find_package(cartesian_interface QUIET)
find_package(catkin REQUIRED COMPONENTS roscpp
                                       std_msgs
                                       geometry_msgs
                                       visualization_msgs
                                       tf)
if(${cartesian_interface_FOUND})
    message("Adding src/CoMStabilizerAddon.cpp to compilation")
    add_library(CoMStabilizerAddon SHARED src/CoMStabilizerAddon.cpp)
    target_link_libraries(CoMStabilizerAddon PRIVATE ${OpenSoT_LIBRARIES}
                                                          CartesianInterface
                                                          ${compliant_stabilizer_LIBRARIES})
    library_install(CoMStabilizerAddon 1 0 0)

endif()

########################################################################
# use YCM to export compliant_stabilizer so that it can be found using find_package #
########################################################################

set(VARS_PREFIX "compliant_stabilizer")

set(compliant_stabilizer_MAJOR_VERSION 0)
set(compliant_stabilizer_MINOR_VERSION 0)
set(compliant_stabilizer_PATCH_VERSION 1)
set(compliant_stabilizer_VERSION ${compliant_stabilizer_MAJOR_VERSION}.${compliant_stabilizer_MINOR_VERSION}.${compliant_stabilizer_PATCH_VERSION})

find_package(YCM REQUIRED)
include(YCMDefaultDirs)
ycm_default_dirs(${VARS_PREFIX})

target_include_directories(compliant_stabilizer PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                                            "$<INSTALL_INTERFACE:${${VARS_PREFIX}_INSTALL_INCLUDEDIR}>")

set_target_properties(compliant_stabilizer PROPERTIES VERSION ${${VARS_PREFIX}_VERSION}
                                SOVERSION ${${VARS_PREFIX}_VERSION})

target_link_libraries(compliant_stabilizer LINK_PUBLIC )

install(DIRECTORY include/ 
        DESTINATION "${${VARS_PREFIX}_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h*")

install(TARGETS compliant_stabilizer  
        EXPORT compliant_stabilizer
        ARCHIVE DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT lib
        RUNTIME DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${${VARS_PREFIX}_INSTALL_LIBDIR}" COMPONENT shlib)
        
#enabling it will add all idynutils dependencies as dependencies for third party users
set_property(GLOBAL APPEND PROPERTY ${VARS_PREFIX}_TARGETS compliant_stabilizer)


include(InstallBasicPackageFiles)
install_basic_package_files(compliant_stabilizer   VARS_PREFIX ${VARS_PREFIX}
                                        VERSION ${${VARS_PREFIX}_VERSION}
                                        COMPATIBILITY SameMajorVersion
                                        TARGETS_PROPERTY ${VARS_PREFIX}_TARGETS
                                        NO_CHECK_REQUIRED_COMPONENTS_MACRO)

include(AddUninstallTarget)
