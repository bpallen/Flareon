# add resources as targets so they show up in IDEs

file(GLOB flareon_shader "${FLAREON_SOURCE_DIR}/res/shader/*")

add_custom_target(
	res
	SOURCES
	${flareon_shader}
)

source_group("shader" FILES ${flareon_shader})
