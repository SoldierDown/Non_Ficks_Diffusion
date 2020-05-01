if(ENABLE_NON_FICKS_DIFFUSION)
  add_executable(converter
    main.cpp
  )

  find_package(OpenMP)
  if(OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
  
  set(OpenVDB_PATH /usr/local/lib/cmake/OpenVDB/)
  list(APPEND CMAKE_MODULE_PATH "${OpenVDB_PATH}")


  find_package(OpenVDB REQUIRED)
  if(OPENVDB_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenVDB_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenVDB_CXX_FLAGS}")
  endif()
  target_link_libraries(converter
    PRIVATE
    nova::dynamics
    nova::geometry
    nova::spgrid
    nova::tools
    OpenVDB::openvdb
    )

  
  file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  
  install(TARGETS converter RUNTIME DESTINATION ${REL})
endif()
