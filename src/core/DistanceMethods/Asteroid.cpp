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


void Asteroid::_getBestSPRRecMissing(unsigned int s,
    StopCriterion stopCriterion,
    std::vector<unsigned int> sprime, 
    std::vector<corax_unode_t *> W0s, 
    corax_unode_t *Wp, 
    corax_unode_t *Wsminus1, 
    corax_unode_t *Vsminus1, 
    std::vector<double> delta_Vsminus2_Wp, // previous deltaAB
    corax_unode_t *Vs, 
    double Lsminus1, // L_s-1
    corax_unode_t *&bestRegraftNode,
    SubBMEToUpdate &subBMEToUpdate,
    double &bestLs,
    unsigned int &bestS,
    const BoolMatrix &belongsToPruned,
    const BoolMatrix &hasChildren,
    std::vector<bool> Vsminus2HasChildren, // does Vsminus2 have children after the previous moves
    std::vector<bool> Vsminus1HasChildren) // does Vsminus1 have children after the previous moves
{
  if (s > stopCriterion.maxRadius || 
      stopCriterion.noImprovement > stopCriterion.maxRadiusWithoutImprovement) {
    return;
  }
  subBMEToUpdate.tempBetween.push_back(Wsminus1);
  corax_unode_t *Ws = getOtherNext(Vs, Vsminus1->back)->back; 
  unsigned int K = _K;
  // compute L_s
  // we compute LS for an NNI to ((A,B),(C,D)) by swapping B and C
  // where:
  // - A is Vsminus1
  // - B is Wp
  // - C is Ws
  // - D is Vs->back
  // A was modified by the previous (simulated) NNI moves
  // and its average distances need an update
  double diff = 0.0;
  std::vector<double> delta_Vsminus1_Wp(K, std::numeric_limits<double>::infinity());
  std::vector<bool> VsHasChildren = Vsminus1HasChildren;
  for (unsigned int k = 0; k < K; ++k) {
    // if Wp has no children in the induced tree, the SPR move has no effect
    // on the score of this family
    if (!hasChildren[Wp->node_index][k]) {
      continue;
    }
    // if there is no regrafting branch in the induced tree anymore, there is no 
    // point in continuing computing stuff for this family
    if (!hasChildren[Vsminus1->back->node_index][k]) {
      continue;
    }
    // The current NNI move affects the current family score if and only if each of
    // its 4 subtrees (after having recursively applied the previous NNI moves!!)
    // have children in the induced tree
    // At this point, Wp has children. Vsminus1HasChildren tells if the updated Vs has children
    // belongsToPruned == true iif both the node children (Ws and Vs->back) have children
    bool applyDiff = Vsminus1HasChildren[k] && belongsToPruned[Vsminus1->back->node_index][k];
 
    // only matters from s == 3
    if (nullptr == W0s[k] && Vsminus1HasChildren[k]) {
      W0s[k] = Wsminus1;
    }

    // deltaCD and deltaBD are trivial because not affected by previous NNI moves
    double deltaCD = getCell(Ws->node_index, Vs->back->node_index, k);
    double deltaBD = getCell(Wp->node_index, Vs->back->node_index, k);
    
    // deltaAB and deltaAC are affected by the previous moves
    double deltaAB = 0.0;
    double deltaAC = 0.0;
    if (sprime[k] <= 1 && applyDiff) {
      deltaAB = getCell(W0s[k]->node_index, Wp->node_index, k);
      deltaAC = getCell(W0s[k]->node_index, Ws->node_index, k);
    } else if (W0s[k]) {
      deltaAC = _pows[sprime[k]] * 
        (getCell(W0s[k]->node_index, Ws->node_index, k) - 
         getCell(Wp->node_index, Ws->node_index, k));
      deltaAC += getCell(Vsminus1->node_index, Ws->node_index, k);
      if (hasChildren[Wsminus1->node_index][k] && Vsminus2HasChildren[k]) {
        deltaAB = 0.5 * (delta_Vsminus2_Wp[k] + 
          getCell(Wsminus1->node_index, Wp->node_index, k));
      } else if (hasChildren[Wsminus1->node_index][k] && !Vsminus2HasChildren[k]) {
        deltaAB = getCell(Wsminus1->node_index, Wp->node_index, k);
      } else if (!hasChildren[Wsminus1->node_index][k] && Vsminus2HasChildren[k]) {
        deltaAB = delta_Vsminus2_Wp[k];
      } else {
        deltaAB = std::numeric_limits<double>::infinity(); // won't be used anyway
      }
    }
    if (applyDiff) {
      double diffk = 0.125 * (deltaAB + deltaCD - deltaAC - deltaBD);
      diff += diffk;
      sprime[k] += 1;
    }
    
    delta_Vsminus1_Wp[k] = deltaAB; // in our out the if???

    // update VsHasChildren for next call
    VsHasChildren[k] = Vsminus1HasChildren[k] || hasChildren[Ws->node_index][k];
  }
  ParallelContext::sumDouble(diff);
  double Ls = Lsminus1 + diff;
  if (Ls > bestLs) {
    bestLs = Ls;
    bestRegraftNode = Vs;
    bestS = s;
    subBMEToUpdate.between = subBMEToUpdate.tempBetween;
    subBMEToUpdate.after = Vs->back;
    subBMEToUpdate.pruned = Wp;
    stopCriterion.noImprovement = 0;
  } else {
    stopCriterion.noImprovement++;
  }
  if (Vs->back->next) {
    _getBestSPRRecMissing(s+1, stopCriterion, sprime, W0s, Wp, Ws, Vs, delta_Vsminus1_Wp, Vs->back->next, 
        Ls, bestRegraftNode, subBMEToUpdate, bestLs, bestS, belongsToPruned, hasChildren, Vsminus1HasChildren, VsHasChildren);
    _getBestSPRRecMissing(s+1, stopCriterion, sprime, W0s, Wp, Ws, Vs, delta_Vsminus1_Wp, Vs->back->next->next, 
        Ls, bestRegraftNode, subBMEToUpdate, bestLs, bestS, belongsToPruned, hasChildren, Vsminus1HasChildren, VsHasChildren);
  }
  subBMEToUpdate.tempBetween.pop_back();
}


void Asteroid::getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves)
{
  _toUpdate.reset();
  for (auto pruneNode: speciesTree.getPostOrderNodes()) {
    double bestDiff = 0.0;
    unsigned int bestS = 1;
    corax_unode_t *bestRegraftNode = nullptr;
    if (getBestSPRFromPrune(maxRadiusWithoutImprovement, pruneNode, bestRegraftNode, bestDiff, bestS)) {
      bestMoves.push_back(SPRMove(pruneNode, bestRegraftNode, bestDiff));
    }
  }
  std::sort(bestMoves.begin(), bestMoves.end());
}



bool Asteroid::getBestSPRFromPrune(unsigned int maxRadiusWithoutImprovement,
      corax_unode_t *prunedNode,
      corax_unode_t *&bestRegraftNode,
      double &bestDiff,
      unsigned int &bestS)
{
  bool foundBetter = false;
  double oldBestDiff = bestDiff;
  auto Wp = prunedNode;
  if (!Wp->back->next) {
    return foundBetter;
  }
  std::vector<corax_unode_t *> V0s;
  V0s.push_back(Wp->back->next->back);
  V0s.push_back(Wp->back->next->next->back);
  unsigned int K = _K;
  std::vector<unsigned int> sprime(K, 1);
  StopCriterion stopCriterion;
  stopCriterion.maxRadiusWithoutImprovement = maxRadiusWithoutImprovement;
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
      std::vector<corax_unode_t *> W0s(K, nullptr);
      std::vector<double> delta_V0_Wp; // not used at first iteration
      std::vector<bool> V0HasChildren = _hasChildren[V0->node_index];
      std::vector<bool> Vminus1HasChildren(V0HasChildren.size()); // won't be read at first iteration
      _getBestSPRRecMissing(s, stopCriterion, sprime, W0s, Wp,  Wsminus1, V,delta_V0_Wp,
          V1, 0.0, bestRegraftNode, _toUpdate, bestDiff, bestS, 
          _belongsToPruned, _hasChildren, Vminus1HasChildren, V0HasChildren);
      if (bestDiff > oldBestDiff) {
        foundBetter = true;
        _toUpdate.before = V0;
        oldBestDiff = bestDiff;
      }
    }
  }
  return foundBetter;
}

