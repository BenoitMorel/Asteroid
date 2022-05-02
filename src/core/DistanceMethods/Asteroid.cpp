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

static DistanceMatrix getNullMatrix(size_t N,
    double value = 0.0)
{
  std::vector<double> nullDistances(N, value);
  return DistanceMatrix(N, nullDistances); 
}

Asteroid::Asteroid(const PLLUnrootedTree &speciesTree, 
      const std::vector<BitVector> &perFamilyCoverage,
      const UIntMatrix &gidToSpid,
      const std::vector<DistanceMatrix> &geneDistanceMatrices):
  _geneDistanceMatrices(geneDistanceMatrices),
  _perFamilyCoverage(perFamilyCoverage)
{
  ParallelContext::barrier();
  _K = geneDistanceMatrices.size();
  _spidToGid.resize(_K);
  _superToInducedNodes.resize(_K);
  _superToInducedNodesRegraft.resize(_K);
  _pruneRegraftDiff.resize(_K);
  for (unsigned int k = 0; k < _K; ++k) {
    _spidToGid[k].resize(speciesTree.getLeavesNumber());
    for (unsigned int gid = 0; gid < gidToSpid[k].size(); ++gid) {
      _spidToGid[k][gidToSpid[k][gid]] = gid;
    }
    _inducedNodeNumber.push_back(gidToSpid[k].size() * 4 - 6);
  }
  for (unsigned int i = 0; i < speciesTree.getLeavesNumber() + 3; ++i) {
    _pows.push_back(std::pow(0.5, i));
  }
  _prunedSpeciesMatrices = std::vector<DistanceMatrix>(_K);
  _avDistances.resize(_K);
  for (unsigned int k = 0; k < _K; ++k) {
    auto geneLeafNumber = gidToSpid[k].size();
    auto geneNodeNumber = _inducedNodeNumber[k];
    _prunedSpeciesMatrices[k] = getNullMatrix(static_cast<size_t>(geneLeafNumber));
    _avDistances[k] = MatrixDouble(geneNodeNumber, 
        std::vector<double>(geneNodeNumber));
    for (unsigned int gid1 = 0; gid1 < geneLeafNumber; ++gid1) {
      for (unsigned int gid2 = 0; gid2 < geneLeafNumber; ++gid2) {
        setCell(gid1, gid2, k, _geneDistanceMatrices[k][gid1][gid2]);
      }
    }
    _pruneRegraftDiff[k] = MatrixDouble(_inducedNodeNumber[k],
        std::vector<double>(_inducedNodeNumber[k], 0.0));
  }
  _computeInducedTrees(speciesTree);
  ParallelContext::barrier();
}


double Asteroid::computeLength(const PLLUnrootedTree &speciesTree)
{
  double res = 0.0;
  _computeInducedTrees(speciesTree);
  for (unsigned k = 0; k < _geneDistanceMatrices.size(); ++k) {
    auto geneLeafNumber = _inducedSpeciesTrees[k]->getLeavesNumber();
    for (unsigned int i = 0; i < geneLeafNumber; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        auto v = _geneDistanceMatrices[k][i][j];
        Logger::info << _prunedSpeciesMatrices[k][i][j] << std::endl;
        res += v * _pows[_prunedSpeciesMatrices[k][i][j]];
      }
    }
  }
  ParallelContext::sumDouble(res);
  return res;
}
  
void Asteroid::_computeInducedTrees(const PLLUnrootedTree &speciesTree)
{
  _inducedSpeciesTrees.clear();
  StringToUint speciesLabelToSpid;
  for (auto leaf: speciesTree.getLeaves()) {
    speciesLabelToSpid.insert({std::string(leaf->label), leaf->node_index});
  }
  for (unsigned int k = 0; k < _K; ++k) {
    _inducedSpeciesTrees.push_back(
        speciesTree.getInducedTree(_perFamilyCoverage[k]));
    StringToUint labelToGid;
    for (auto leaf: _inducedSpeciesTrees[k]->getLeaves()) {
      auto spid = speciesLabelToSpid[leaf->label];
      auto gid = _spidToGid[k][spid];
      labelToGid.insert({std::string(leaf->label), gid});
    }
    _inducedSpeciesTrees[k]->reindexLeaves(labelToGid);
  }
  auto superPostOrderNodes = speciesTree.getPostOrderNodes();
  for (unsigned int k = 0; k < _K; ++k) {
    speciesTree.mapNodesWithInducedTree(*_inducedSpeciesTrees[k],
        superPostOrderNodes,
        _superToInducedNodes[k],
        _superToInducedNodesRegraft[k]);
    InternodeDistance::computeFromSpeciesTree(*_inducedSpeciesTrees[k],
         _prunedSpeciesMatrices[k]);
  }
}

void Asteroid::_computeAvDistances()
{
  for (unsigned int k = 0; k < _K; ++k) {
    auto &prunedSpeciesTree = *_inducedSpeciesTrees[k];
    auto subtrees1 = prunedSpeciesTree.getReverseDepthNodes();
    for (auto n1: subtrees1) {
      auto i1 = n1->node_index;
      auto subtrees2 = prunedSpeciesTree.getPostOrderNodesFrom(n1->back);
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
    double diffMinus1) // L_s-1
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
    deltaAC -= _pows[s] * getCell(Wp->node_index, Ws->node_index, k);
    deltaAC += _pows[s] * getCell(W0->node_index, Ws->node_index, k);
  }
  double diff = diffMinus1 + 0.125 * (deltaAB + deltaCD - deltaAC - deltaBD);
  _pruneRegraftDiff[k][Wp->node_index][Vs->node_index] = diff;
  _pruneRegraftDiff[k][Wp->node_index][Vs->back->node_index] = diff;
  // recursive call
  if (Vs->back->next) {
    precomputeSPRDiffRec(k, s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next, 
        diff);
    precomputeSPRDiffRec(k, s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next->next, 
        diff);
  }
}

void Asteroid::precomputeSPRDiffFromPrune(unsigned int k, 
    corax_unode_t *pruneNode)
{
  auto Wp = pruneNode;
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
          V1, 0.0);
    }
  }
}



void Asteroid::getBestSPRFromPruneRec(StopCriterion stopCriterion,
    corax_unode_t *regraft,
    double lastBestScore,
    const std::vector<unsigned int> &ks,
    const std::vector<unsigned int> &inducedPruneIndices,
    corax_unode_t *&bestRegraft,
    double &bestScore)
{
  // stop if we went too far without improvement
  if (stopCriterion.noImprovement > 
      stopCriterion.maxRadiusWithoutImprovement) {
    return;
  }
  double newScore = 0.0;
  // iterate over all induced species tree on which the pruned
  // subtree is non empty
  for (unsigned int i = 0; i < ks.size(); ++i) {
    auto k = ks[i];
    auto inducedPruneIndex = inducedPruneIndices[i];
    auto iRegraft = _superToInducedNodesRegraft[k][regraft->node_index]; 
    // diff is the diff of the score after applying this SPR
    // move to this induced species tree
    auto diff =  _pruneRegraftDiff[k][inducedPruneIndex][iRegraft->node_index];
    newScore += diff;  
  }
  ParallelContext::sumDouble(newScore);
  
  if (newScore > lastBestScore) {
    // Improvement, search further
    stopCriterion.noImprovement = 0;
    lastBestScore = newScore;
  } else {
    stopCriterion.noImprovement++;
  }
  if (newScore > bestScore) {
    bestScore = newScore;
    bestRegraft = regraft;
  }
  if (regraft->next) {
    auto left = PLLUnrootedTree::getLeft(regraft);
    auto right = PLLUnrootedTree::getRight(regraft);
    getBestSPRFromPruneRec(stopCriterion, left, lastBestScore, ks,inducedPruneIndices,  bestRegraft, bestScore);
    getBestSPRFromPruneRec(stopCriterion, right, lastBestScore, ks, inducedPruneIndices,  bestRegraft, bestScore);
  }
}

bool Asteroid::getBestSPRFromPrune(unsigned int maxRadiusWithoutImprovement,
    corax_unode_t *pruneNode,
    corax_unode_t *&bestRegraftNode,
    double &bestDiff)
{
  bool foundBetter = false;
  double oldBestDiff = bestDiff;
  
  if (!pruneNode->back->next) {
    return foundBetter;
  }
  
  std::vector<unsigned int> ks;
  std::vector<unsigned int> inducedPruneIndices;
  for (unsigned int k = 0; k < _K; ++k) {
    auto inducedPrune = _superToInducedNodes[k][pruneNode->node_index];
    if (inducedPrune) {
      ks.push_back(k);
      inducedPruneIndices.push_back(inducedPrune->node_index);
    }
  }

  std::vector<corax_unode_t *> regrafts;
  auto r1 = pruneNode->back->next->back;
  auto r2 = pruneNode->back->next->next->back;
  if (r1->next) {
    regrafts.push_back(r1);
  }
  if (r2->next) {
    regrafts.push_back(r2);
  }
  for (auto regraft: regrafts) {
    StopCriterion stopCriterion;
    stopCriterion.maxRadiusWithoutImprovement = maxRadiusWithoutImprovement;
    getBestSPRFromPruneRec(stopCriterion, 
        regraft,
        0.0,
        ks,
        inducedPruneIndices,
        bestRegraftNode,
        bestDiff);
    if (bestDiff > oldBestDiff) {
      foundBetter = true;
      oldBestDiff = bestDiff;
    }
  }
  return foundBetter;
}

    
void Asteroid::getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves)
{
  _computeAvDistances();
  bestMoves.clear();
  for (unsigned int k = 0; k < _K; ++k) {
    for (auto pruneNode: _inducedSpeciesTrees[k]->getPostOrderNodes()) {
      auto &regraftDiff = _pruneRegraftDiff[k][pruneNode->node_index];
      std::fill(regraftDiff.begin(), regraftDiff.end(), 0.0);
      precomputeSPRDiffFromPrune(k, pruneNode);
    }
  }
  for (auto pruneNode: speciesTree.getPostOrderNodes()) {
    double bestDiff = 0.0;
    corax_unode_t *bestRegraftNode = nullptr;
    if (getBestSPRFromPrune(maxRadiusWithoutImprovement, 
          pruneNode, 
          bestRegraftNode, 
          bestDiff)) {
      bestMoves.push_back(SPRMove(pruneNode, bestRegraftNode, bestDiff));
    }
  }
  std::sort(bestMoves.begin(), bestMoves.end());
}



