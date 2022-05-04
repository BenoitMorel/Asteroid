#pragma once
#include <map>
#include <memory>
#include <util/types.hpp>
#include <DistanceMethods/AsteroidTypes.hpp>
#include <trees/InducedStepwiseTree.hpp>
#include <trees/StepwiseTree.hpp>
#include <trees/InducedStepwiseTree.hpp>
#include <set>

using InducedStepwiseTrees = std::vector<std::unique_ptr<InducedStepwiseTree> >;

using IntPairToDouble = std::map <std::pair<unsigned int, unsigned int>, double>;

class StepwiseAsteroid {
public:
  StepwiseAsteroid(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells);

  void insertLabel(const std::string &label);

  void exportTree(const std::string &outputPath);
  
private:
  size_t _N; // total number of taxa
  size_t _K; // number of gene trees (for this rank)
  const StringToUint &_speciesToSpid;
  std::vector<std::string> _spidToSpecies;
  const std::vector<GeneCell *> &_geneCells;
  UIntMatrix _spidToGid;
  StepwiseTree _tree;
  InducedStepwiseTrees _inducedTrees;
  BitVector _insertedSpecies;
  std::vector<MatrixUint> _speciesMatrices;
  std::vector<double> _pows;

  // information to compute the score diff
  std::vector<MatrixDouble> _deltas;
  std::vector<IntPairToDouble> _slowDeltas;
  MatrixDouble _phis;
  
  void _updateDeltas();
  void _updatePhis(unsigned int spid);
  double _computeDelta(std::set<Node *>&s1,
      std::set<Node *>s2,
      const DistanceMatrix &geneMatrix,
      const MatrixUint &speciesMatrix,
      const std::vector<unsigned int> &spidToGid);
  double _getPhi(unsigned int gid,
      Node *node,
      const DistanceMatrix &geneMatrix,
      const std::vector<unsigned int> &spidToGid,
      unsigned int depth);
};



