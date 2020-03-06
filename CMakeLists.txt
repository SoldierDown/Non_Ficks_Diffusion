set(ENABLE_NON_FICKS_DIFFUSION OFF CACHE BOOL "Enable building the Non-Fick's Diffusion project tree.")

cmake_policy(SET CMP0054 NEW)

if(ENABLE_NON_FICKS_DIFFUSION)

add_subdirectory(nfd_2d)
add_subdirectory(nfd_3d)
add_subdirectory(Ficks_Solver)
# add_subdirectory(plugin_Non_Ficks_Diffusion)
add_subdirectory(plugin_Heightfield)
endif(ENABLE_NON_FICKS_DIFFUSION)
