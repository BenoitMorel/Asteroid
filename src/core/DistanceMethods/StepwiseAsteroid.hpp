#pragma once
#include <map>
#include <memory>
#include <util/types.hpp>
#include <DistanceMethods/AsteroidTypes.hpp>
#include <trees/InducedStepwiseTree.hpp>
#include <trees/StepwiseTree.hpp>
#include <trees/InducedStepwiseTree.hpp>

using InducedStepwiseTrees = std::vector<std::unique_ptr<InducedStepwiseTree> >;

using IntPairToDouble = std::map <std::pair<unsigned int, unsigned int>, double>;

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
  UIntMatrix _spidToGid;
  StepwiseTree _tree;
  InducedStepwiseTrees _inducedTrees;
  BitVector _insertedSpecies;
  std::vector<DistanceMatrix> _speciesMatrices;

  // information to compute the score diff
  std::vector<IntPairToDouble> _deltas;
  MatrixDouble _phis;

  void _updateDeltas(unsigned int spid);
};



