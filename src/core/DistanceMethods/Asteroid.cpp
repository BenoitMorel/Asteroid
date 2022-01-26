#include "Asteroid.hpp"
#include <IO/Logger.hpp>
#include <limits>
#include <algorithm>
#include <DistanceMethods/InternodeDistance.hpp>

using DistanceVectorMatrix = std::vector<DistanceMatrix>;

static corax_unode_t *getOtherNext(corax_unode_t *n1, 
    corax_unode_t *n2)
{
  if (n1->next == n2) {
    return n2->next;
  } else {
    return n1->next;
  }
}

static BoolMatrix getBoolMatrix(unsigned int N,
    unsigned int M,
    bool value = false)
{
  std::vector<bool> v(M, value);
  return BoolMatrix(N, v);
}

static DistanceMatrix getNullMatrix(unsigned int N,
    double value = 0.0)
{
  std::vector<double> nullDistances(N, value);
  return DistanceMatrix(N, nullDistances); 
}
static DistanceMatrix getMatrix(unsigned int N,
    unsigned int M,
    double value)
{
  std::vector<double> nullDistances(M, value);
  return DistanceMatrix(N, nullDistances); 
}


/**
 *  Computes the distances between a leaf (corresponding to
 *  currentNode at the first recursive call) and all other 
 *  nodes in the unrooted tree. If belongsToPruned is set,
 *  the distances are computed on the pruned subtree
 */
static void fillDistancesRec(corax_unode_t *currentNode, 
    double currentDistance,
    std::vector<double> &distances,
    std::vector<bool> *belongsToPruned)
{
  if (!belongsToPruned || (*belongsToPruned)[currentNode->node_index]) {
    currentDistance += 1.0;
  }
  if (!currentNode->next) {
    // leaf
    if (!belongsToPruned || (*belongsToPruned)[currentNode->node_index]) {
      auto spid = reinterpret_cast<intptr_t>(currentNode->data);
      distances[spid] = currentDistance;
    }
    return;
  }
  fillDistancesRec(currentNode->next->back, 
      currentDistance, 
      distances,
      belongsToPruned);
  fillDistancesRec(currentNode->next->next->back, 
      currentDistance, 
      distances,
      belongsToPruned);
}

/**
 *  Fills the internode distance matrix of a pruned species tree.
 *  belongsToPruned tells which species leaves do not belong
 *  to the pruned species tree.
 */
static void fillSpeciesDistances(const PLLUnrootedTree &speciesTree,
    DistanceMatrix &speciesDistanceMatrix,
    std::vector<bool> *belongsToPruned = nullptr)
{
  for (auto leaf: speciesTree.getLeaves()) {
    if (belongsToPruned && !(*belongsToPruned)[leaf->node_index]) {
      continue;
    }
    std::string label = leaf->label;
    auto id = reinterpret_cast<intptr_t>(leaf->data);
    fillDistancesRec(leaf->back,
        0.0,
        speciesDistanceMatrix[id],
        belongsToPruned);
  }
}

static void getPrunedSpeciesMatrix(const PLLUnrootedTree &speciesTree,
  std::vector<bool> &coverage,
  DistanceMatrix &distanceMatrix)
{
  std::vector<bool> hasChildren(speciesTree.getDirectedNodesNumber(), false);
  std::vector<bool> belongsToPruned(speciesTree.getDirectedNodesNumber(), false);
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto index = node->node_index;
    if (node->next) {
      auto leftIndex = node->next->back->node_index;
      auto rightIndex = node->next->next->back->node_index;
      hasChildren[index] = hasChildren[leftIndex] || hasChildren[rightIndex];
      belongsToPruned[index] = hasChildren[leftIndex] && hasChildren[rightIndex];
    } else {
      auto spid = reinterpret_cast<intptr_t>(node->data);
      bool isCovered = coverage[spid];
      hasChildren[index] = isCovered;
      belongsToPruned[index] = isCovered;
    }
  }
  fillSpeciesDistances(speciesTree,
      distanceMatrix,
      &belongsToPruned);
}
    


Asteroid::Asteroid(const PLLUnrootedTree &speciesTree, 
      const BoolMatrix &perFamilyCoverage,
      const UIntMatrix &gidToSpid,
      const std::vector<DistanceMatrix> &geneDistanceMatrices):
  _gidToSpid(gidToSpid),
  _geneDistanceMatrices(geneDistanceMatrices),
  _perFamilyCoverage(perFamilyCoverage)
{
  ParallelContext::barrier();
  _K = geneDistanceMatrices.size();
  _N = speciesTree.getDirectedNodesNumber(); 
  _speciesNumber = speciesTree.getLeavesNumber();
  _spidToGid.resize(_K);
  _superToInducedNodes.resize(_K);
  _inducedToSuperNodes.resize(_K);
  _pruneRegraftDiff.resize(_K);
  for (unsigned int k = 0; k < _K; ++k) {
    _spidToGid[k].resize(_N);
    for (unsigned int gid = 0; gid < _gidToSpid[k].size(); ++gid) {
      _spidToGid[k][_gidToSpid[k][gid]] = gid;
    }
  }
  for (unsigned int i = 0; i < _N + 3; ++i) {
    _pows.push_back(std::pow(0.5, i));
  }
  _prunedSpeciesMatrices = std::vector<DistanceMatrix>(_K);
  _subBMEs.resize(_K);
  for (unsigned int k = 0; k < _K; ++k) {
    auto geneNumbers = gidToSpid[k].size();
    auto directedGeneNumbers = geneNumbers * 4 - 6;
    _prunedSpeciesMatrices[k] = getNullMatrix(geneNumbers);
    _subBMEs[k] = std::vector<double>(directedGeneNumbers * directedGeneNumbers, 
      std::numeric_limits<double>::infinity());
    for (unsigned int i = 0; i < geneNumbers; ++i) {
      for (unsigned int j = 0; j < geneNumbers; ++j) {
        setCell(i, j, k, _geneDistanceMatrices[k][i][j]);
      }
    }
    _pruneRegraftDiff[k] = MatrixDouble(directedGeneNumbers,
        std::vector<double>(directedGeneNumbers, 0.0));
  }
  ParallelContext::barrier();
}


double Asteroid::computeBME(const PLLUnrootedTree &speciesTree)
{
  double res = 0.0;
  _prunedSpeciesTrees.clear();
  StringToUint speciesLabelToSpid;
  for (auto leaf: speciesTree.getLeaves()) {
    speciesLabelToSpid.insert({std::string(leaf->label), leaf->node_index});
  }
   
  for (unsigned int k = 0; k < _K; ++k) {
    _prunedSpeciesTrees.push_back(
        speciesTree.getInducedTree(_perFamilyCoverage[k]));
    StringToUint labelToGid;
    for (auto leaf: _prunedSpeciesTrees[k]->getLeaves()) {
      auto spid = speciesLabelToSpid[leaf->label];
      auto gid = _spidToGid[k][spid];
      labelToGid.insert({std::string(leaf->label), gid});
    }
    _prunedSpeciesTrees[k]->reindexLeaves(labelToGid);
    speciesTree.mapNodesWithInducedTree(*_prunedSpeciesTrees[k],
        _superToInducedNodes[k],
        _inducedToSuperNodes[k]);
    InternodeDistance::computeFromSpeciesTree(*_prunedSpeciesTrees[k],
         _prunedSpeciesMatrices[k]);
  }
  for (unsigned k = 0; k < _geneDistanceMatrices.size(); ++k) {
    auto geneNumbers = _gidToSpid[k].size();
    for (unsigned int i = 0; i < geneNumbers; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        auto v = _geneDistanceMatrices[k][i][j];
        res += v * _pows[_prunedSpeciesMatrices[k][i][j]];
      }
    }
  }
  ParallelContext::sumDouble(res);
  _computeSubBMEsPrune(speciesTree);
  Logger::timed << "after compute pruned " << std::endl;
  return res;
}



void Asteroid::_computeSubBMEsPrune(const PLLUnrootedTree &speciesTree)
{
  for (unsigned int k = 0; k < _K; ++k) {
    auto &prunedSpeciesTree = *_prunedSpeciesTrees[k];
    auto subtrees1 = prunedSpeciesTree.getReverseDepthNodes();
    auto nodesNumber = subtrees1.size();
    for (auto n1: subtrees1) {
      auto subtrees2 = prunedSpeciesTree.getPostOrderNodesFrom(n1->back);
      auto i1 = n1->node_index;
      for (auto n2: subtrees2) {
        auto i2 = n2->node_index;
        if (!n1->next && !n2->next) { // both leaves
          // already done
        } else if (n1->next && !n2->next) { // n2 is a leaf
          // we already computed the symetric
          setCell(i1, i2, k, getCell(i2, i1, k));
        } else  { // n2 is not a leaf
          auto left2 = n2->next->back->node_index;
          auto right2 = n2->next->next->back->node_index;
          auto v = 0.5 * (getCell(i1, left2, k) + 
              getCell(i1, right2, k));
          setCell(i1, i2, k, v);
        }
      }
    }
  }
}


// test regrafting  Wp between Vs and Vsplus1 after 
// having regrafted Wp between Vsminus and Vs
// and recursively test further SPR mvoes
void Asteroid::precomputeSPRDiffRec(unsigned int k,
    unsigned int s,
    corax_unode_t *W0, 
    corax_unode_t *Wp, 
    corax_unode_t *Wsminus1, 
    corax_unode_t *Vsminus1, 
    double delta_Vsminus2_Wp, // previous deltaAB
    corax_unode_t *Vs, 
    double diffMinus1, // L_s-1
    std::vector<double> &regraftDiff)
{
  corax_unode_t *Ws = getOtherNext(Vs, Vsminus1->back)->back; 
  // compute L_s
  // we compute LS for an NNI to ((A,B),(C,D)) by swapping B and C
  // where:
  // - A is Vsminus1
  // - B is Wp
  // - C is Ws
  // - D is Vs->back
  // A was modified by the previous (simulated) NNI moves
  // and its average distances need an update
  double deltaCD = getCell(Ws->node_index, Vs->back->node_index, k);
  double deltaBD = getCell(Wp->node_index, Vs->back->node_index, k);
  double deltaAB = 0.0;
  double deltaAC = 0.0;
  if (s == 1) {
    deltaAB = getCell(W0->node_index, Wp->node_index, k);
    deltaAC = getCell(W0->node_index, Ws->node_index, k);
  } else {
    deltaAB = 0.5 * (delta_Vsminus2_Wp + 
      getCell(Wsminus1->node_index, Wp->node_index, k));
    deltaAC = getCell(Vsminus1->node_index, Ws->node_index, k);
    deltaAC -= pow(0.5, s) * getCell(Wp->node_index, Ws->node_index, k);
    deltaAC += pow(0.5, s) * getCell(W0->node_index, Ws->node_index, k);
  }
  double diff = diffMinus1 + 0.125 * (deltaAB + deltaCD - deltaAC - deltaBD);
  //Logger::info << regraftDiff.size() << " " << Vs->node_index << std::endl;
  regraftDiff[Vs->node_index] = diff;
  // recursive call
  if (Vs->back->next) {
    precomputeSPRDiffRec(k, s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next, 
        diffMinus1, regraftDiff);
    precomputeSPRDiffRec(k, s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next->next, 
        diffMinus1, regraftDiff);
  }
}

void Asteroid::precomputeSPRDiffFromPrune(unsigned int k, 
    corax_unode_t *prunedNode,
    std::vector<double> &regraftDiff)
{
  auto Wp = prunedNode;
  if (!Wp->back->next) {
    return;
  }
  std::vector<corax_unode_t *> V0s;
  V0s.push_back(Wp->back->next->back);
  V0s.push_back(Wp->back->next->next->back);
  for (auto V0: V0s) {
    auto V = getOtherNext(Wp->back, V0->back);
    if (!V->back->next) {
      continue;
    }
    std::vector<corax_unode_t *> V1s;
    V1s.push_back(V->back->next);
    V1s.push_back(V->back->next->next);
    for (auto V1: V1s) {
      unsigned int s = 1;
      corax_unode_t *W0 = V0;
      corax_unode_t *Wsminus1 = W0;
      double deltaAB = 0.0; // not used at first iteration
      precomputeSPRDiffRec(k, s, W0, Wp,  Wsminus1, V, deltaAB,
          V1, 0.0, regraftDiff);
    }
  }
}

void Asteroid::getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves)
{
  _toUpdate.reset();
  for (unsigned int k = 0; k < _K; ++k) {
    for (auto pruneNode: _prunedSpeciesTrees[k]->getPostOrderNodes()) {
        precomputeSPRDiffFromPrune(k, pruneNode, _pruneRegraftDiff[k][pruneNode->node_index]);
    }
  }
  Logger::info <<"getBestSPR " << std::endl;
  unsigned int i = 0;
  for (auto prunedNode: speciesTree.getPostOrderNodes()) {
    double bestPrunedDiff = 0.0;
    auto p = prunedNode->node_index;
    corax_unode_t *bestRegraftNode = nullptr;
    for (auto regraftNode: speciesTree.getPostOrderNodesFrom(prunedNode->back)) {
      auto r = regraftNode->node_index;
      double diff = 0.0;
      for (unsigned int k = 0; k < _K; ++k) {
        auto inducedPrune = _superToInducedNodes[k][p];
        auto inducedRegraft = _superToInducedNodes[k][r];
        if (inducedPrune && inducedRegraft) {
          diff += _pruneRegraftDiff[k][inducedPrune->node_index][inducedRegraft->node_index];
        }
      }
      ParallelContext::sumDouble(diff);
      if (diff > bestPrunedDiff) {
        bestPrunedDiff = diff;
        bestRegraftNode = regraftNode;
      }
    }
    if (bestRegraftNode) {
      bestMoves.push_back(SPRMove(prunedNode, bestRegraftNode, bestPrunedDiff));
    }
  }
  std::sort(bestMoves.begin(), bestMoves.end());
  /*
  Logger::info << "Number of bestMoves: " << bestMoves.size() << std::endl;
  for (auto &m: bestMoves) {
    Logger::info << m.score << std::endl;
  }
  */

}




