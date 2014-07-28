
# allow targets to include from each other
include_directories("${FLAREON_SOURCE_DIR}/src")

# this needs to be done first
# insert common lib here

add_subdirectory("./flareon")
