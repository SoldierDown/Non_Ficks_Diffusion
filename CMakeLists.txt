set(ENABLE_NON_FICKS_DIFFUSION OFF CACHE BOOL "Enable building the Non-Fick's Diffusion project tree.")

cmake_policy(SET CMP0054 NEW)

if(ENABLE_NON_FICKS_DIFFUSION)

add_subdirectory(mpm_2d)
add_subdirectory(mpm_3d)
add_subdirectory(plugin_Non_Ficks_Diffusion)

endif(ENABLE_NON_FICKS_DIFFUSION)
