#pragma once

#include <memory>
#include <vector>
#include <DistanceMethods/Asteroid.hpp>


class AsteroidOptimizer {
public:
  AsteroidOptimizer(PLLUnrootedTree &speciesTree,
      const BoolMatrix &perFamilyCoverage,
      const UIntMatrix &gidToSpid,
      const std::vector<DistanceMatrix> &distanceMatrices);
  virtual ~AsteroidOptimizer() {}
  double optimize();
private:
  double eval(PLLUnrootedTree &tree);
  bool computeAndApplyBestSPR();

  PLLUnrootedTree &_speciesTree;
  Asteroid _asteroid;
  double _lastScore;
};



