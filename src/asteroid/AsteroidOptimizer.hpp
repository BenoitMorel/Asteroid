#pragma once

#include <memory>
#include <vector>
#include <DistanceMethods/Asteroid.hpp>

/**
 *  Optimize a species tree from a set of gene distance matrices 
 *  and their mappings to the species, using the Asteroid
 *  algorithm
 */
class AsteroidOptimizer {
public:
  /**
   *  Init the optimizer
   */
  AsteroidOptimizer(PLLUnrootedTree &speciesTree,
      const BoolMatrix &perFamilyCoverage,
      const UIntMatrix &gidToSpid,
      const std::vector<DistanceMatrix> &distanceMatrices);
  virtual ~AsteroidOptimizer() {}
  /**
   *  Runs the optimization
   */
  double optimize();
private:
  /**
   * Evaluates the score of the current tree
   */
  double eval(PLLUnrootedTree &tree);
  
  /**
   * Evaluates candidate SPR moves, and 
   * applies several of them simultaneously
   * Returns true if the score improved
   */
  bool computeAndApplyBestSPR();

  PLLUnrootedTree &_speciesTree;
  Asteroid _asteroid;
  double _lastScore;
};



