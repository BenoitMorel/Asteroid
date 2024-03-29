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

  // information to compute the score diff

  // _deltas[k][i1][i2] stores the average distance between
  // the nodes whose indices are i1 and i2 in the induced tree k
  // See the FastME paper for the definition of the average distance
  std::vector<MatrixDouble> _deltas;
  
  void _updateDeltas();
  double _computeDelta(const std::set<Node *>&s1,
      const std::set<Node *> &s2,
      const DistanceMatrix &geneMatrix,
      const MatrixUint &speciesMatrix,
      const std::vector<unsigned int> &spidToGid);
  double _getPhi(unsigned int gid,
      Node *node,
      const DistanceMatrix &geneMatrix,
      const std::vector<unsigned int> &spidToGid,
      unsigned int depth);
  double _updateDeltasNewTaxonAux(
    unsigned int taxonGid,
    unsigned int taxonIndex,
    unsigned int k,
    Node *node);
  void _updateDeltasNewTaxon(unsigned int taxonSpid);
  void _updateDeltasAux(unsigned int k,
      Node *insertedNode,
      Node *B,
      std::vector<Node *> &chain);
  void _updateDeltasChain(unsigned int k,
    Node *newTaxon,
    const std::vector<Node *> &chain);
};



