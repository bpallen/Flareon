# add resources as targets so they show up in IDEs

file(GLOB flareon_shader "${FLAREON_SOURCE_DIR}/res/shader/*")
file(GLOB flareon_lensdata "${FLAREON_SOURCE_DIR}/res/lensdata/*")

add_custom_target(
	res
	SOURCES
	${flareon_shader}
	${flareon_lensdata}
)

source_group("shader" FILES ${flareon_shader})
source_group("lensdata" FILES ${flareon_lensdata})
