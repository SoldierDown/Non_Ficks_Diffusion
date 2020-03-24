set(ENABLE_NON_FICKS_DIFFUSION_PLUGIN OFF CACHE BOOL "Enable non-Fick's diffusion plugin.")
if(ENABLE_NON_FICKS_DIFFUSION_PLUGIN)
  find_package(assimp)

  if(NOT assimp_FOUND)
    message(FATAL_ERROR "Cannot build Non-Fick's Diffusion plugin without the Assimp library being present.")
  endif()
  
  add_library(plugin_Non_Ficks_Diffusion
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Model.cpp
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Mesh.cpp
    main.cpp
    )

  if(WIN32)
    target_compile_definitions(plugin_Non_Ficks_Diffusion PUBLIC OPENGL_3D_WIN32=1)
  endif()
  
  if(UNIX)
    target_compile_definitions(plugin_Non_Ficks_Diffusion PUBLIC OPENGL_3D_LINUX=1)
  endif()

  target_include_directories(plugin_Non_Ficks_Diffusion
                             PUBLIC ${VISUALIZER_KERNEL_INCLUDES})
  
  
  target_link_libraries(plugin_Non_Ficks_Diffusion
    PRIVATE
    nova::dynamics
    nova::geometry
    nova::spgrid
    nova::tools
    ${Boost_LIBRARIES}
    ${assimp_LIBRARIES}
    plugin_Assimp_miniball
    opengl_kernel
    )
  
  target_compile_features(plugin_Non_Ficks_Diffusion
    PRIVATE
    cxx_auto_type
    )
endif()
