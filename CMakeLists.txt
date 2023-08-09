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

cmake_minimum_required(VERSION 3.0)
PROJECT(compliant_stabilizer)

find_package(catkin)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(OpenSoT QUIET)

INCLUDE_DIRECTORIES(include ${EIGEN3_INCLUDE_DIR})

set(compliant_stabilizer_sources src/compliantstabilizer.cpp
                          src/FilterHClass.cpp)

if (${OpenSoT_FOUND})
    message("Adding src/CoMStabilizer.cpp to compilation")
    set(compliant_stabilizer_sources ${compliant_stabilizer_sources}
        src/CoMStabilizer.cpp)

    INCLUDE_DIRECTORIES(${OpenSoT_INCLUDE_DIRS})

endif()

catkin_package(

    INCLUDE_DIRS

    LIBRARIES
    compliant_stabilizer

    )


ADD_LIBRARY(compliant_stabilizer SHARED ${compliant_stabilizer_sources})

if (${OpenSoT_FOUND})
    TARGET_LINK_LIBRARIES(compliant_stabilizer ${OpenSoT_LIBRARIES})
endif()

install(TARGETS compliant_stabilizer
        DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
        )

install(DIRECTORY include DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
FILES_MATCHING PATTERN "*.h")









