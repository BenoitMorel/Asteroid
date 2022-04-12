#pragma once

#include <util/types.hpp>
#include <DistanceMethods/AsteroidTypes.hpp>
#include <trees/StepwiseTree.hpp>


class StepwiseAsteroid {
public:
  StepwiseAsteroid(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells);

  void insertLabel(const std::string &label);

  void exportTree(const std::string &outputPath);
  
private:
  unsigned int _N; // total number of taxa
  unsigned int _K; // number of gene trees (for this rank)
  const StringToUint &_speciesToSpid;
  std::vector<std::string> _spidToSpecies;
  const std::vector<GeneCell *> &_geneCells;
  StepwiseTree _tree;
  const std::vector<StepwiseTree> _inducedTrees;
  BitVector _insertedSpecies;
  DistanceMatrix _speciesMatrix;
};



