set(ENABLE_HEIGHTFIELD_PLUGIN OFF CACHE BOOL "Enable heightfield plugin.")

if(ENABLE_HEIGHTFIELD_PLUGIN)
  find_package(assimp)

  if(NOT assimp_FOUND)
    message(FATAL_ERROR "Cannot build embedded deformables plugin without the Assimp library being present.")
  endif()
  
  add_library(plugin_Heightfield
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Model.cpp
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Mesh.cpp
    main.cpp
    )

  if(WIN32)
    target_compile_definitions(plugin_Heightfield PUBLIC OPENGL_3D_WIN32=1)
  endif()
  
  if(UNIX)
    target_compile_definitions(plugin_Heightfield PUBLIC OPENGL_3D_LINUX=1)
  endif()

  target_include_directories(plugin_Heightfield
                             PUBLIC ${VISUALIZER_KERNEL_INCLUDES})
  
  
  target_link_libraries(plugin_Heightfield
    PRIVATE
    nova::geometry
    nova::tools
    ${Boost_LIBRARIES}
    ${assimp_LIBRARIES}
    plugin_Assimp_miniball
    opengl_kernel
    )
  
  target_compile_features(plugin_Heightfield
    PRIVATE
    cxx_auto_type
    )
endif()
