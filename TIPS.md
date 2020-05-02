# CMakelists
find_package(TBB)
  find_package(Blosc)
  find_package(IlmBase)
  find_package(OpenVDB)
  
  target_link_libraries(converter
    PUBLIC
    openvdb
    Half
    tbb
  )
# dependencies
- uninstall libilmbase and reinstall it

# zip
-  find [0-9]  [0-3][0-9] | zip -r part.zip -@