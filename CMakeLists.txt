project(openvdb_visualizer)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_definitions(-std=c++14)
if(WIN32)
    add_definitions(-W4)
    add_definitions(-DNOMINMAX=1)
else(WIN32)
    add_definitions(-Wall)
    add_definitions(-wx)
endif(WIN32)

#add_subdirectory(extension)
add_subdirectory(plugin)
#add_subdirectory(shader)

install(FILES README.md
	    DESTINATION docs/)
