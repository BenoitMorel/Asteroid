
set(asteroid-core_SOURCES 
  IO/Logger.cpp
  maths/Random.cpp
  parallelization/ParallelContext.cpp
  )

add_library(asteroid-core STATIC ${asteroid-core_SOURCES})

target_include_directories(asteroid-core
  PUBLIC ${ASTEROID_INCLUDE_DIRS}
  )


