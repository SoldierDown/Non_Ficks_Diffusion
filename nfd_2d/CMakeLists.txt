if(ENABLE_NON_FICKS_DIFFUSION)
  add_executable(nfd_2d
    main.cpp
    ../MPM_Driver.cpp
    ../MPM_Example.cpp
  )

  find_package(OpenMP)
  if(OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()

  target_link_libraries(nfd_2d
    PRIVATE
    nova::dynamics
    nova::geometry
    nova::spgrid
    nova::tools
    )
  
  file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  
  install(TARGETS nfd_2d RUNTIME DESTINATION ${REL})
endif()
