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

/**
 *  Implements the global prune length evaluation for a given
 *  species tree and its SPR neighbors, according to the
 *  Asteroid algorithm
 */
class Asteroid {
public:
  Asteroid(const PLLUnrootedTree &speciesTree, 
      const BoolMatrix &perFamilyCoverage,
      const UIntMatrix &gidToSpid,
      const std::vector<DistanceMatrix> &geneDistanceMatrices);
  virtual ~Asteroid() {}
  /*
   *  Computes the BME score of the species tree
   *  and update the _avDistances used to speedup the
   *  precomputeSPRDiff call.
   */
  virtual double computeLength(const PLLUnrootedTree &speciesTree);

  /**
   *  Tries all possible SPR moves from the speciesTree
   *  and returns the best one as well as the difference
   *  between the current score and the new score
   */
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves);


private:
  struct StopCriterion {
    unsigned int maxRadius;
    unsigned int maxRadiusWithoutImprovement;
    unsigned int noImprovement;
    StopCriterion(): maxRadius(99999),
      maxRadiusWithoutImprovement(4),
      noImprovement(0)
      {}
  };
  // Index convention:
  // - k -> gene tree index
  // - gid: node_index in the induced species tree
  // - spid: node_index in the species tree

  // _gidToSpid[k][gid] == spid
  const UIntMatrix &_gidToSpid;
  // _spidToGid[k][spid] == gid
  UIntMatrix _spidToGid;
  // internode gene distance matrices
  const std::vector<DistanceMatrix> &_geneDistanceMatrices;
  // number of families 
  size_t _K;
  // _perFamilyCoverage[k][i] is true if the family k covers the species i
  BoolMatrix _perFamilyCoverage;
  std::vector<unsigned int> _inducedNodeNumber;
  // _pows[i] == pow(2, i) (precomputed to speedup computations)
  std::vector<double> _pows;
  // trees induced by the current species tree and
  // each family coverage pattern
  std::vector<std::shared_ptr<PLLUnrootedTree> > _inducedSpeciesTrees; 
  // internode distance matrices of the induced species trees
  std::vector<DistanceMatrix> _prunedSpeciesMatrices;
  // mapping between species tree prining nodes and 
  // induced tree pruning nodes
  std::vector<NodeVector> _superToInducedNodes; 
  // mapping between species tree regrafting nodes and 
  // induced tree regrafting nodes
  std::vector<NodeVector> _superToInducedNodesRegraft;
  // average distance between every pair of species subtrees 
  // Precomputing it allows to evaluate SPR moves in
  // constant time
  // this is the Delta matrix in the FastME2 supp mat paper
  //
  std::vector<MatrixDouble> _avDistances;
  // _pruneRegraftDiff[k][i][j] stores the score diff
  // resulting in pruning the node with index i and regrafting
  // it at the node j in the induced tree k
  std::vector<MatrixDouble> _pruneRegraftDiff;
  
  
  
  double getCell(size_t sp1, size_t sp2, size_t k) {
    return _avDistances[k][sp1][sp2];
  }
  void setCell(size_t sp1, size_t sp2, size_t k, double v) {
    _avDistances[k][sp1][sp2] = v;
  }

  void getBestSPRFromPruneRec(StopCriterion stopCriterion,
    corax_unode_t *regraft,
    double lastScore,
    const std::vector<unsigned int> &ks,
    const std::vector<unsigned int> &inducedPruneIndices,
    corax_unode_t *&bestRegraft,
    double &bestScore);
  



  void _computeInducedTrees(const PLLUnrootedTree &speciesTree);
  void _computeAvDistances();

  void precomputeSPRDiffFromPrune(unsigned int k, 
      corax_unode_t *prunedNode);

  void precomputeSPRDiffRec(unsigned int k,
      unsigned int s,
    corax_unode_t *W0, 
    corax_unode_t *Wp, 
    corax_unode_t *Wsminus1, 
    corax_unode_t *Vsminus1, 
    double delta_Vsminus2_Wp, // previous deltaAB
    corax_unode_t *Vs, 
    double diffMinus1); // L_s-1

  bool getBestSPRFromPrune(unsigned int maxRadiusWithoutImprovement,
    corax_unode_t *pruneNode,
    corax_unode_t *&bestRegraftNode,
    double &bestDiff);
};
