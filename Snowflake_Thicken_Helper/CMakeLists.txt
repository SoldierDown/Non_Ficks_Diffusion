if(ENABLE_NON_FICKS_DIFFUSION)
  add_executable(thicken_helper
    main.cpp
  )

  find_package(OpenMP)
  if(OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()

  find_package(TBB)
  find_package(Blosc)
  find_package(IlmBase)
  find_package(OpenVDB)
  
  target_link_libraries(thicken_helper
    PUBLIC
    openvdb
    Half
    tbb
  )
	target_link_libraries(thicken_helper
    PRIVATE
    nova::dynamics
    nova::geometry
    nova::spgrid
    nova::tools
    )
  
  	file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  
  	install(TARGETS converter RUNTIME DESTINATION ${REL})
endif()
