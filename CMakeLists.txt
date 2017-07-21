project(openvdb_visualizer)
cmake_minimum_required(VERSION 2.8)

option(BUILD_MAYA_PLUGIN "Build the maya plugin." ON)
option(BUILD_MTOA_EXTENSION "Build the MtoA extension." OFF)
option(BUILD_ARNOLD_SHADER "Build the Arnold shader." OFF)
option(INCLUDE_HEADERS_IN_BUILD "Include the headers next to the source files for IDEs that need this." OFF)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (WIN32)
    add_definitions(-W4)
    add_definitions(-DNOMINMAX=1)
else ()
    add_definitions(-Wall)
    add_definitions(-wx)
endif ()

add_definitions(-std=c++11)

if (BUILD_MAYA_PLUGIN)
	add_subdirectory(plugin)
endif ()
if (BUILD_MTOA_EXTENSION)
	add_subdirectory(extension)
endif ()
if (BUILD_ARNOLD_SHADER)
	add_subdirectory(shader)
endif ()

install(FILES README.md
	    DESTINATION docs/)
