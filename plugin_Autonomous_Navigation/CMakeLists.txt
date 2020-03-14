if(ENABLE_AUTONOMOUS_NAVIGATION)
  find_package(assimp)

  if(NOT assimp_FOUND)
    message(FATAL_ERROR "Cannot build embedded deformables plugin without the Assimp library being present.")
  endif()

  add_library(plugin_Autonomous_Navigation
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Model.cpp
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Mesh.cpp
    main.cpp
    )

  if(WIN32)
    target_compile_definitions(plugin_Autonomous_Navigation PUBLIC OPENGL_3D_WIN32=1)
  endif()
  
  if(UNIX)
    target_compile_definitions(plugin_Autonomous_Navigation PUBLIC OPENGL_3D_LINUX=1)
  endif()

  target_include_directories(plugin_Autonomous_Navigation
                             PUBLIC ${VISUALIZER_KERNEL_INCLUDES})
  
  
  target_link_libraries(plugin_Autonomous_Navigation
    PRIVATE
    nova::dynamics
    nova::spgrid
    nova::tools
    ${Boost_LIBRARIES}
    ${assimp_LIBRARIES}
    plugin_Assimp_miniball
    )
  
  target_compile_features(plugin_Autonomous_Navigation
    PRIVATE
    cxx_auto_type
    )
endif()
