#pragma once

#include <util/types.hpp>
#include <DistanceMethods/AsteroidTypes.hpp>
#include <trees/StepwiseTree.hpp>


class StepwiseAsteroid {
public:
  StepwiseAsteroid(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells,
    const std::array<std::string, 3> initialLabels);

  void insertLabel(const std::string &label);

  void exportTree(const std::string &outputPath);
  
private:
  unsigned int _N; // total number of taxa
  unsigned int _K; // number of gene trees (for this rank)
  const StringToUint &_speciesToSpid;
  std::vector<std::string> _spidToSpecies;
  const std::vector<GeneCell *> &_geneCells;
  StepwiseTree _tree;
  BitVector _insertedSpecies;
  DistanceMatrix _speciesMatrix;
};



