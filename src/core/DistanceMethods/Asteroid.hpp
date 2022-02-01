#pragma once

#include <string>
#include <util/types.hpp>
#include <vector>
#include <unordered_set>
#include <trees/PLLUnrootedTree.hpp>


using BoolMatrix = std::vector< std::vector<bool> >;
using UIntMatrix = std::vector< std::vector<unsigned int> >;

struct SPRMove {
  corax_unode_t *pruneNode;
  corax_unode_t *regraftNode;
  double score;
  SPRMove(corax_unode_t *pruneNode,
      corax_unode_t *regraftNode,
      double score): 
    pruneNode(pruneNode),
    regraftNode(regraftNode),
    score(score) {}
  
  bool operator < (const SPRMove& other) const {
    return other.score < score;
  }
  bool operator ==(const SPRMove& other) const {
    return other.score == score;
  }
};

class Asteroid {
public:
  Asteroid(const PLLUnrootedTree &speciesTree, 
      const BoolMatrix &perFamilyCoverage,
      const UIntMatrix &gidToSpid,
      const std::vector<DistanceMatrix> &geneDistanceMatrices);
  virtual ~Asteroid() {}
  /*
   *  Computes the BME score of the species tree
   *  and update the _subBMEs used to speedup the
   *  precomputeSPRDiff call.
   */
  virtual double computeBME(const PLLUnrootedTree &speciesTree);

  /**
   *  Tries all possible SPR moves from the speciesTree
   *  and returns the best one as well as the difference
   *  between the current score and the new score
   */
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves);

  void _computedInducedTrees(const PLLUnrootedTree &speciesTree);
  void _updateInducedTrees(const PLLUnrootedTree &speciesTree);

private:
  // Internal implementation
  //
  // Index convention:
  // - k -> gene tree index
  // - gid: node_index in the induced species tree
  // - spid: node_index in the species tree
  //
  //

  // _gidToSpid[k][gid] == spid
  const UIntMatrix &_gidToSpid;
  // _spidToGid[k][spid] == gid
  UIntMatrix _spidToGid;
  // _precomputeSPRDiffRecMissing[k][i][j] is the distance
  // between species i and j for the gene family k
  // this value is only valid if no inf 
  const std::vector<DistanceMatrix> &_geneDistanceMatrices;
  // number of families 
  size_t _K;
  // _perFamilyCoverage[k][i] is true if the family k covers the species i
  BoolMatrix _perFamilyCoverage;
  std::vector<unsigned int> _inducedNodeNumber;

  // _prunedSpeciesMatrices[k] is the internode distance for the
  // species tree induced by the family k
  std::vector<std::shared_ptr<PLLUnrootedTree> > _inducedSpeciesTrees; 
  std::vector<NodeVector> _superToInducedNodes; 
  std::vector< std::vector<NodeSet> > _inducedToSuperNodes; 
  std::vector< std::vector<NodeSet> > _inducedToSuperNodesRegraft; 
  std::vector<DistanceMatrix> _prunedSpeciesMatrices;
  std::vector<MatrixDouble> _pruneRegraftDiff;
  std::vector<MatrixDouble> _subBMEs;
  double getCell(size_t sp1, size_t sp2, size_t k) {
    return _subBMEs[k][sp1][sp2];
  }
  void setCell(size_t sp1, size_t sp2, size_t k, double v) {
    _subBMEs[k][sp1][sp2] = v;
  }
  
  // _hasChildren[i][k] == true if there is at least one leaf
  // under i that belongs to the species tree induced by the
  // family k
  BoolMatrix _hasChildren;
  // _belongsToPruned[i][k] == true if the node i belongs to
  // the species tree induced by the family k
  BoolMatrix _belongsToPruned;
  // _pows[i] == pow(2, i) (precomputed to speedup computations)
  std::vector<double> _pows;

  double _computeBMEPrune(const PLLUnrootedTree &speciesTree);
  
  void _computeSubBMEsPrune();

  struct StopCriterion {
    unsigned int maxRadius;
    unsigned int maxRadiusWithoutImprovement;
    unsigned int noImprovement;
    StopCriterion(): maxRadius(99999),
      maxRadiusWithoutImprovement(4),
      noImprovement(0)
      {}
  };

  void precomputeSPRDiffFromPrune(unsigned int k, 
      corax_unode_t *prunedNode,
    std::vector<double> &regraftDiff);

  void precomputeSPRDiffRec(unsigned int k,
      unsigned int s,
    corax_unode_t *W0, 
    corax_unode_t *Wp, 
    corax_unode_t *Wsminus1, 
    corax_unode_t *Vsminus1, 
    double delta_Vsminus2_Wp, // previous deltaAB
    corax_unode_t *Vs, 
    double diffMinus1, // L_s-1
    std::vector<double> &regraftDiff);
};
