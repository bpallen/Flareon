
# get source files
# we could list these manually...
file(GLOB flareon_src "*.cpp" "*.c")
file(GLOB flareon_hdr "*.hpp" "*.h")

add_executable(flareon ${flareon_src} ${flareon_hdr})

set_target_properties(
	flareon
    PROPERTIES
    LINKER_LANGUAGE CXX
)

# TODO export glew definitions properly
add_definitions(-DGLEW_STATIC -DGLEW_MX)
target_link_libraries(flareon glew glfw ${GLFW_LIBRARIES})
