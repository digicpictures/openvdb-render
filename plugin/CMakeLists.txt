set(PLUGIN_NAME openvdb_render)

file(GLOB SRC *.cpp)
file(GLOB HDR_h *.h)
file(GLOB HDR_hpp *.hpp)

configure_file(shader/volume.ogsfx ${CMAKE_BINARY_DIR}/bin/shader/volume.ogsfx COPYONLY)

add_library(${PLUGIN_NAME} SHARED ${SRC} ${HDR_h} ${HDR_hpp})
set_target_properties(${PLUGIN_NAME} PROPERTIES PREFIX "")
if (WIN32)
    set_target_properties(${PLUGIN_NAME} PROPERTIES SUFFIX ".mll")
endif ()
target_link_libraries(${PLUGIN_NAME} ${CONAN_LIBS})

install(TARGETS ${PLUGIN_NAME}
        DESTINATION plug-ins)

# cgfx shader for Maya 2015, ogsfx for Maya 2017.
install(FILES shader/volume.cgfx DESTINATION plug-ins/shader)
install(FILES shader/volume.ogsfx DESTINATION plug-ins/shader)

file(GLOB TEMPLATES template/*.py)

install(FILES ${TEMPLATES} DESTINATION scripts)
