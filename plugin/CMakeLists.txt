set(PLUGIN_NAME openvdb_render)

# Compile-time reduction
include(U:/development/gilianz/common/cmake/cotire.cmake)

if(WIN32)
    add_definitions(-D_BOOL -DREQUIRE_IOSTREAM)
endif(WIN32)

# Maya
find_package(Maya 2015)
include_directories(${MAYA_INCLUDE_DIRS})
link_directories(${MAYA_LIBRARY_DIR})
set(MAYA_LIBS OpenMaya Foundation OpenMayaAnim OpenMayaFX OpenMayaUI OpenMayaRender)

# OpenGL
find_package(OpenGL REQUIRED)
add_library(GL INTERFACE)
set_target_properties(GL   PROPERTIES INTERFACE_LINK_LIBRARIES      ${OPENGL_gl_LIBRARY})
add_library(GLU INTERFACE)
set_target_properties(GLU  PROPERTIES INTERFACE_LINK_LIBRARIES      ${OPENGL_glu_LIBRARY})

if(WIN32)
    set(DIGIC_DEV_LIB_DIR U:/development/_shared/libs)
else(WIN32)
    message(FATAL_ERROR "TODO: set DIGIC_DEV_LIB_DIR for linux.")
endif(WIN32)

# GLEW
add_library(GLEW SHARED IMPORTED)
set(GLEW_ROOT ${DIGIC_DEV_LIB_DIR}/glew/2.0.0)
set_target_properties(GLEW PROPERTIES IMPORTED_LOCATION_RELEASE     ${GLEW_ROOT}/lib_win64/glew32.dll)
set_target_properties(GLEW PROPERTIES IMPORTED_IMPLIB_RELEASE       ${GLEW_ROOT}/lib_win64/glew32.lib)
set_target_properties(GLEW PROPERTIES IMPORTED_LOCATION_DEBUG       ${GLEW_ROOT}/lib_win64/glew32d.dll)
set_target_properties(GLEW PROPERTIES IMPORTED_IMPLIB_DEBUG         ${GLEW_ROOT}/lib_win64/glew32d.lib)
set_target_properties(GLEW PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${GLEW_ROOT}/include)

# Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.59 REQUIRED COMPONENTS thread)
get_filename_component(BOOST_ROOT ${Boost_INCLUDE_DIRS} DIRECTORY)
link_directories(${BOOST_ROOT}/lib_win64)

# TBB
add_library(TBB SHARED IMPORTED)
set(TBB_ROOT ${DIGIC_DEV_LIB_DIR}/TBB/4.0)
set_target_properties(TBB PROPERTIES IMPORTED_LOCATION_RELEASE     ${TBB_ROOT}/bin_win64/tbb.dll)
set_target_properties(TBB PROPERTIES IMPORTED_IMPLIB_RELEASE       ${TBB_ROOT}/lib_win64/tbb.lib)
set_target_properties(TBB PROPERTIES IMPORTED_LOCATION_DEBUG       ${TBB_ROOT}/bin_win64/tbb_debug.dll)
set_target_properties(TBB PROPERTIES IMPORTED_IMPLIB_DEBUG         ${TBB_ROOT}/lib_win64/tbb_debug.lib)
set_target_properties(TBB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${TBB_ROOT}/include)
set_target_properties(TBB PROPERTIES INTERFACE_LINK_LIBRARIES Boost::thread)

# OpenEXR::Half
add_library(OpenEXR_Half STATIC IMPORTED)
set(OpenEXR_ROOT ${DIGIC_DEV_LIB_DIR}/OpenEXR/2.2.0_static)
set_target_properties(OpenEXR_Half PROPERTIES IMPORTED_LOCATION_RELEASE     ${OpenEXR_ROOT}/lib_win64/half_static_vc140_md.lib)
set_target_properties(OpenEXR_Half PROPERTIES IMPORTED_LOCATION_DEBUG       ${OpenEXR_ROOT}/lib_win64/half_static_vc140_md_gd.lib)
set_target_properties(OpenEXR_Half PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${OpenEXR_ROOT}/include)

# ZLIB
add_library(ZLIB STATIC IMPORTED)
set(ZLIB_ROOT ${DIGIC_DEV_LIB_DIR}/zlib/1.2.8)
set_target_properties(ZLIB PROPERTIES IMPORTED_LOCATION_RELEASE     ${ZLIB_ROOT}/lib_win64/zlibstatic.lib)
set_target_properties(ZLIB PROPERTIES IMPORTED_LOCATION_DEBUG       ${ZLIB_ROOT}/lib_win64/zlibstatic.lib)
set_target_properties(ZLIB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${ZLIB_ROOT}/include)

# Blosc
add_library(Blosc SHARED IMPORTED)
set(Blosc_ROOT ${DIGIC_DEV_LIB_DIR}/blosc)
set_target_properties(Blosc PROPERTIES IMPORTED_LOCATION_RELEASE     ${Blosc_ROOT}/lib_win64/blosc.dll)
set_target_properties(Blosc PROPERTIES IMPORTED_IMPLIB_RELEASE       ${Blosc_ROOT}/lib_win64/blosc.lib)
set_target_properties(Blosc PROPERTIES IMPORTED_LOCATION_DEBUG       ${Blosc_ROOT}/lib_win64/blosc.dll)
set_target_properties(Blosc PROPERTIES IMPORTED_IMPLIB_DEBUG         ${Blosc_ROOT}/lib_win64/blosc.lib)
set_target_properties(Blosc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${Blosc_ROOT}/include)

# OpenVDB
add_library(OpenVDB STATIC IMPORTED)
set(OPENVDB_ROOT ${DIGIC_DEV_LIB_DIR}/OpenVDB/3.2.0)
set_target_properties(OpenVDB PROPERTIES IMPORTED_LOCATION_RELEASE      ${OPENVDB_ROOT}/lib_win64/OpenVDB_static_vc140_md.lib)
set_target_properties(OpenVDB PROPERTIES IMPORTED_LOCATION_DEBUG        ${OPENVDB_ROOT}/lib_win64/OpenVDB_static_vc140_md_gd.lib)
set_target_properties(OpenVDB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${OPENVDB_ROOT}/include;${OPENVDB_ROOT}/include/openvdb")
set_target_properties(OpenVDB PROPERTIES INTERFACE_COMPILE_DEFINITIONS "OPENVDB_STATICLIB;OPENVDB_OPENEXR_STATICLIB;OPENEXR_USE_BLOSC")
set_target_properties(OpenVDB PROPERTIES INTERFACE_LINK_LIBRARIES      "Blosc;OpenEXR_Half;TBB;ZLIB")

# lib
file(GLOB SRC *.cpp)
file(GLOB HDR_h *.h)
file(GLOB HDR_hpp *.hpp)

add_library(${PLUGIN_NAME} SHARED ${SRC} ${HDR_h} ${HDR_hpp} volume.cgfx)
set_target_properties(${PLUGIN_NAME} PROPERTIES PREFIX "")
if (WIN32)
    set_target_properties(${PLUGIN_NAME} PROPERTIES SUFFIX ".mll")
endif ()
target_link_libraries(${PLUGIN_NAME} ${MAYA_LIBS} GLEW GL OpenVDB)

cotire(${PLUGIN_NAME})

install(TARGETS ${PLUGIN_NAME}
        DESTINATION plug-ins)
install(FILES volume.cgfx
        DESTINATION plug-ins/shader/)

file(GLOB TEMPLATES template/*.py)

install(FILES ${TEMPLATES} DESTINATION scripts)
