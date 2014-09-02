file(GLOB gfft_src "*.cpp")
file(GLOB gfft_inc "${GFFT_SOURCE_DIR}/include/GFFT/*.h")

add_library(gfft ${gfft_src} ${gfft_inc})


set_target_properties(gfft
    PROPERTIES
    LINKER_LANGUAGE CXX
) 
