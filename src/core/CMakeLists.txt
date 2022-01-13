
set(asteroid-core_SOURCES 
  DistanceMethods/InternodeDistance.cpp
  DistanceMethods/Asteroid.cpp
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  maths/Random.cpp
  trees/PLLUnrootedTree.cpp
  parallelization/ParallelContext.cpp
  )

add_library(asteroid-core STATIC ${asteroid-core_SOURCES})

target_include_directories(asteroid-core
  PUBLIC ${ASTEROID_INCLUDE_DIRS}
  )


