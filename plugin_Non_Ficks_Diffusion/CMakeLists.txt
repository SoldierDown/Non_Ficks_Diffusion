if(ENABLE_NON_FICKS_DIFFUSION)
  find_package(assimp)

  if(NOT assimp_FOUND)
    message(FATAL_ERROR "Cannot build Non-Fick's Diffusion plugin without the Assimp library being present.")
  endif()
  
  add_library(plugin_MPM
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Model.cpp
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Mesh.cpp
    main.cpp
    )

  if(WIN32)
    target_compile_definitions(plugin_MPM PUBLIC OPENGL_3D_WIN32=1)
  endif()
  
  if(UNIX)
    target_compile_definitions(plugin_MPM PUBLIC OPENGL_3D_LINUX=1)
  endif()

  target_include_directories(plugin_MPM
                             PUBLIC ${VISUALIZER_KERNEL_INCLUDES})
  
  
  target_link_libraries(plugin_MPM
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
  
  target_compile_features(plugin_MPM
    PRIVATE
    cxx_auto_type
    )
endif()
