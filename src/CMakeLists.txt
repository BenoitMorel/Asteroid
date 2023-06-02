
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")


add_subdirectory(core)
  
macro(add_program program_name program_sources)
  #if(LINUX) 
  #  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  #  set(BUILD_SHARED_LIBS OFF)
  #  set(CMAKE_EXE_LINKER_FLAGS "-static")
  #endif()
  add_executable(${program_name} ${program_sources})
  target_link_libraries(${program_name} 
    asteroid-core
    ${CORAX_LIBRARIES}
    ${MPI_LIBRARIES}
    )
  set_target_properties(${program_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    )
  target_include_directories(${program_name}
    PUBLIC ${ASTEROID_INCLUDE_DIRS}
    )
endmacro()

add_subdirectory(asteroid)
add_subdirectory(test)


