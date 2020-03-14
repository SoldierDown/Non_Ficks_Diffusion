add_executable(poisson_solver
  main.cpp
  )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(poisson_solver
  PRIVATE
  nova::dynamics
  nova::tools
  nova::geometry
  nova::spgrid
  )

file( RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

install(TARGETS poisson_solver RUNTIME DESTINATION ${REL})
