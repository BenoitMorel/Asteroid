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


static DistanceMatrix getNullMatrix(unsigned int N,
    double value = 0.0)
{
  std::vector<double> nullDistances(N, value);
  return DistanceMatrix(N, nullDistances); 
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
    _inducedDirectedNodeNumbers.push_back(_gidToSpid[k].size() * 4 - 6);
  }
  for (unsigned int i = 0; i < _N + 3; ++i) {
    _pows.push_back(std::pow(0.5, i));
  }
  _prunedSpeciesMatrices = std::vector<DistanceMatrix>(_K);
  _subBMEs.resize(_K);
  for (unsigned int k = 0; k < _K; ++k) {
    auto geneNumbers = _gidToSpid[k].size();
    _prunedSpeciesMatrices[k] = getNullMatrix(geneNumbers);
    auto subBMESize = _inducedDirectedNodeNumbers[k] * 
      _inducedDirectedNodeNumbers[k];
    _subBMEs[k] = std::vector<double>(subBMESize, 
      std::numeric_limits<double>::infinity());
    for (unsigned int i = 0; i < geneNumbers; ++i) {
      for (unsigned int j = 0; j < geneNumbers; ++j) {
        setCell(i, j, k, _geneDistanceMatrices[k][i][j]);
      }
    }
    _pruneRegraftDiff[k] = MatrixDouble(_inducedDirectedNodeNumbers[k],
        std::vector<double>(_inducedDirectedNodeNumbers[k], 0.0));
  }
  ParallelContext::barrier();
}


double Asteroid::computeBME(const PLLUnrootedTree &speciesTree)
{
  double res = 0.0;
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
    speciesTree.mapNodesWithInducedTree(*_inducedSpeciesTrees[k],
        _superToInducedNodes[k],
        _inducedToSuperNodes[k]);
    InternodeDistance::computeFromSpeciesTree(*_inducedSpeciesTrees[k],
         _prunedSpeciesMatrices[k]);
  }
  for (unsigned k = 0; k < _geneDistanceMatrices.size(); ++k) {
    auto geneNumbers = _inducedSpeciesTrees[k]->getLeavesNumber();
    for (unsigned int i = 0; i < geneNumbers; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        auto v = _geneDistanceMatrices[k][i][j];
        res += v * _pows[_prunedSpeciesMatrices[k][i][j]];
      }
    }
  }
  ParallelContext::sumDouble(res);
  _computeSubBMEsPrune();
  Logger::timed << "after compute pruned " << std::endl;
  return res;
}



void Asteroid::_computeSubBMEsPrune()
{
  for (unsigned int k = 0; k < _K; ++k) {
    auto &prunedSpeciesTree = *_inducedSpeciesTrees[k];
    auto subtrees1 = prunedSpeciesTree.getReverseDepthNodes();
    for (auto n1: subtrees1) {
      auto subtrees2 = prunedSpeciesTree.getPostOrderNodesFrom(n1->back);
      auto i1 = n1->node_index;
      for (auto n2: subtrees2) {
        auto i2 = n2->node_index;
        if (!n1->next && !n2->next) { // both leaves
          // already done
        } else if (n1->next && !n2->next) { // n2 is a leaf
          // we already computed the symetric
          assert(! (k == 0 && i2 == 1 && i1 == 0));
          assert(! (k == 0 && i1 == 1 && i2 == 0));
          setCell(i1, i2, k, getCell(i2, i1, k));
        } else  { // n2 is not a leaf
          auto left2 = n2->next->back->node_index;
          auto right2 = n2->next->next->back->node_index;
          auto v = 0.5 * (getCell(i1, left2, k) + 
              getCell(i1, right2, k));
          setCell(i1, i2, k, v);
          assert(! (k == 0 && i2 == 1 && i1 == 0));
          assert(! (k == 0 && i1 == 1 && i2 == 0));
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
        diff, regraftDiff);
    precomputeSPRDiffRec(k, s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next->next, 
        diff, regraftDiff);
  }
}

void Asteroid::precomputeSPRDiffFromPrune(unsigned int k, 
    corax_unode_t *prunedNode,
    std::vector<double> &regraftDiff)
{
  std::fill(regraftDiff.begin(), regraftDiff.end(), 0.0);
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
  bestMoves.clear();
  for (unsigned int k = 0; k < _K; ++k) {
    for (auto pruneNode: _inducedSpeciesTrees[k]->getPostOrderNodes()) {
      precomputeSPRDiffFromPrune(k, pruneNode, _pruneRegraftDiff[k][pruneNode->node_index]);
    }
  }
  auto postOrderNodes = speciesTree.getPostOrderNodes();
  MatrixDouble totalDiff = MatrixDouble(postOrderNodes.size(),
      std::vector<double>(postOrderNodes.size(), 0.0));

  for (unsigned int k = 0; k < _K; ++k) {
    auto inducedPostOrderNodes = 
      _inducedSpeciesTrees[k]->getPostOrderNodes();
    for (auto inducedPrunedNode: inducedPostOrderNodes) {
      auto ipIndex = inducedPrunedNode->node_index;
      auto &superPrunedNodes = _inducedToSuperNodes[k][ipIndex];
      for (auto inducedRegraftNode: inducedPostOrderNodes) {
        auto irIndex = inducedRegraftNode->node_index;
        auto &superRegraftNodes = _inducedToSuperNodes[k][irIndex];
        auto inducedDiff = _pruneRegraftDiff[k][ipIndex][irIndex];
        if (inducedDiff != 0.0) {
          for (auto superPrunedNode: superPrunedNodes) {
            auto spIndex = superPrunedNode->node_index;
            for (auto superRegraftNode: superRegraftNodes) {
              auto srIndex = superRegraftNode->node_index;
              totalDiff[spIndex][srIndex] += inducedDiff;
            }
          }
        }
      }
    }
  }
  for (auto pruneNode: postOrderNodes) {
    auto p = pruneNode->node_index;
    double bestPruneDiff = 0.0;
    corax_unode_t *bestRegraftNode = nullptr;
    for (auto regraftNode: postOrderNodes) {
      auto r = regraftNode->node_index;
      ParallelContext::sumDouble(totalDiff[p][r]);
      if (totalDiff[p][r] > bestPruneDiff) {
        bestPruneDiff = totalDiff[p][r];
        bestRegraftNode = regraftNode;
      }
    }
    if (bestRegraftNode) {
      bestMoves.push_back(SPRMove(pruneNode, bestRegraftNode, bestPruneDiff));
    }
  }
  std::sort(bestMoves.begin(), bestMoves.end());
  /*
   * Logger::info << "Number of bestMoves: " << bestMoves.size() << std::endl;
  for (auto &m: bestMoves) {
    Logger::info << m.score << std::endl;
  }
*/

}




