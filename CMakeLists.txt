project(openvdb_visualizer)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#add_definitions(-Wall)
#add_definitions(-wx) # -Werror for MSVC
add_definitions(-std=c++14)
add_definitions(-DNOMINMAX=1)

#add_subdirectory(extension)
add_subdirectory(plugin)
#add_subdirectory(shader)

install(FILES README.md
	    DESTINATION docs/)
