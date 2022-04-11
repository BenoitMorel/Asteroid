
set(asteroid-core_SOURCES 
  DistanceMethods/InternodeDistance.cpp
  DistanceMethods/Asteroid.cpp
  DistanceMethods/AsteroidOptimizer.cpp
  DistanceMethods/StepwiseAsteroid.cpp
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  IO/ParallelOfstream.cpp
  maths/bitvector.cpp
  maths/Random.cpp
  trees/PLLUnrootedTree.cpp
  trees/SplitHashtable.cpp
  trees/StepwiseTree.cpp
  parallelization/ParallelContext.cpp
  )

add_library(asteroid-core STATIC ${asteroid-core_SOURCES})

target_include_directories(asteroid-core
  PUBLIC ${ASTEROID_INCLUDE_DIRS}
  )


