#
# Copyright (C) 2016, Patrick Seemann
#
# This software may be modified and distributed under the terms
# of the BSD 3-Clause license. See the LICENSE.txt file for details.

cmake_minimum_required(VERSION 3.0)
project(curvest)

include(ExternalProject)

set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED yes)

# Enable OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Optimization flags
IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")
ELSE()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

# Add MVE
find_package(JPEG REQUIRED)
find_package(TIFF REQUIRED)
find_package(PNG REQUIRED)
ExternalProject_Add(MVE
    GIT_REPOSITORY https://github.com/simonfuhrmann/mve.git
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND "")

# Include MVE
ExternalProject_Get_Property(MVE SOURCE_DIR)
set(MVE_LIBRARIES ${SOURCE_DIR}/libs/mve/libmve.a
${SOURCE_DIR}/libs/util/libmve_util.a)
include_directories(${SOURCE_DIR}/libs)

# Add colormaps
ExternalProject_Add(colormap
    GIT_REPOSITORY https://github.com/seepa/colormap.git
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
ExternalProject_Get_Property(colormap SOURCE_DIR)
include_directories(${SOURCE_DIR}/libcolormap)

# add source files
file(GLOB SOURCE_FILES
    "${PROJECT_SOURCE_DIR}/libcurvest/*.cc"
    "${PROJECT_SOURCE_DIR}/libcurvest/*.h"
)
# create static curvest lib
add_library(curvest STATIC ${SOURCE_FILES})
include_directories("${PROJECT_SOURCE_DIR}")

# tools: curvestimate
add_executable(curvestimate "${PROJECT_SOURCE_DIR}/tools/curvestimate.cc")
add_dependencies(curvestimate MVE)
target_link_libraries(curvestimate curvest ${MVE_LIBRARIES} ${PNG_LIBRARIES}
${JPEG_LIBRARIES} ${TIFF_LIBRARIES})

# tools: densityfield
add_executable(densityfield "${PROJECT_SOURCE_DIR}/tools/densityfield.cc")
add_dependencies(densityfield MVE)
target_link_libraries(densityfield curvest ${MVE_LIBRARIES} ${PNG_LIBRARIES}
${JPEG_LIBRARIES} ${TIFF_LIBRARIES})

# tools: mesh colorize
add_executable(mesh_colorize "${PROJECT_SOURCE_DIR}/tools/mesh_colorize.cc")
add_dependencies(mesh_colorize MVE colormap)
target_link_libraries(mesh_colorize ${MVE_LIBRARIES} ${PNG_LIBRARIES}
${JPEG_LIBRARIES} ${TIFF_LIBRARIES})
