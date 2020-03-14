add_executable(smoke_3d
  main.cpp
  ../Smoke_Driver.cpp
  ../Smoke_Example.cpp
  )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(smoke_3d
  PRIVATE
  nova::dynamics
  nova::geometry
  nova::spgrid
  nova::tools
  )

file( RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

install(TARGETS smoke_3d RUNTIME DESTINATION ${REL})
