#include "StepwiseAsteroid.hpp"

#include <limits>
#include <IO/ParallelOfstream.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

StepwiseAsteroid::StepwiseAsteroid(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells):
  _N(speciesToSpid.size()),
  _K(geneCells.size()),
  _speciesToSpid(speciesToSpid),
  _spidToSpecies(_N),
  _geneCells(geneCells),
  _spidToGid(_K),
  _insertedSpecies(_N, false),
  _speciesMatrices(_K),
  _deltas(_K),
  _slowDeltas(_K),
  _phis(_K)

{
  for (auto pair: speciesToSpid) {
    auto label = pair.first;
    auto spid = pair.second;
    _spidToSpecies[spid] = label;
  }
  for (unsigned int k = 0; k < _K; ++k) {
    _inducedTrees.push_back(std::make_unique<InducedStepwiseTree>(
          _tree,
          geneCells[k]->coverage));
    _tree.addListener(_inducedTrees.back().get());

    auto Nk = geneCells[k]->coverage.count();
    auto internalNodes = Nk * 4 - 6;
    _phis[k].resize(internalNodes);

    _deltas[k] = MatrixDouble(internalNodes, 
        std::vector<double>(internalNodes, 0.0));
    _speciesMatrices[k] = UIntMatrix(Nk, 
        std::vector<unsigned int>(Nk, 0));
    auto &spidToGid = _spidToGid[k];
    auto &gidToSpid = _geneCells[k]->gidToSpid;
    spidToGid.resize(_N);
    for (unsigned int gid = 0; gid < gidToSpid.size(); ++gid) {
      spidToGid[gidToSpid[gid]] = gid;
    }
  }
  for (unsigned int i = 0; i < _N; ++i) {
    _pows.push_back(std::pow(0.5, double(i)));
  }
}
  


double getAsteroidScore(const MatrixUint &speciesMatrix,
    const GeneCell *geneCell,
    const BitVector &addedSpids) // the spids that belong to the tree being built
{
  double score = 0.0;
  if ((addedSpids & geneCell->coverage).count() < 4) {
    return score;
  }
  const auto &geneMatrix = geneCell->distanceMatrix;
  const auto &gidToSpid = geneCell->gidToSpid;
  for (unsigned int gid1 = 0; gid1 < geneMatrix.size(); ++gid1) {
    auto spid1 = gidToSpid[gid1];
    if (!addedSpids[spid1]) {
      continue; 
    }
    for (unsigned int gid2 = 0;  gid2 < gid1; ++gid2) {
      auto spid2 = gidToSpid[gid2];
      if (!addedSpids[spid2]) {
        continue; 
      }
      double d = geneMatrix[gid1][gid2];
      auto m = speciesMatrix[gid1][gid2];
      score += d * pow(0.5, m);
    }
  }
  return score / 2.0;
}

/**
 *  Compute the internode distance bitween spid and the leaves under node
 *  Also fill subgids with the gid of the leaves under node
 */
void updateDistanceMatrixAux(unsigned int gid,
    Node *node,
    unsigned int distance,
    std::vector<unsigned int> &spidToGid,
    MatrixUint &distanceMatrix,
    std::set<unsigned int> &subgids)
{
  if (!node->next) {
    auto gid2 = spidToGid[node->spid];
    distanceMatrix[gid][gid2] = distance;
    distanceMatrix[gid2][gid] = distance;
    subgids.insert(gid2);
    return;
  }
  auto left = node->next->back;
  auto right = node->next->next->back;
  updateDistanceMatrixAux(gid, left, distance + 1, spidToGid, distanceMatrix, subgids);
  updateDistanceMatrixAux(gid, right, distance + 1, spidToGid, distanceMatrix, subgids);
}
    

void printMatrix(const DistanceMatrix &m) {
  for (const auto &v: m) {
    for (const auto val: v) {
      Logger::info << val << " ";
    }
    Logger::info << std::endl;
  }
}

static MatrixUint getUpdatedDistanceMatrix(unsigned int gid,
    Node *insertionBranch,
    std::vector<unsigned int> &spidToGid,
    const MatrixUint &speciesMatrix
    )
{
  MatrixUint res(speciesMatrix);
  if (!insertionBranch || !insertionBranch->back) {
    // edge case:
    // this is the second inserted taxa
    // we don't need to update the matrix
    return res;
  }
  auto leftSubtree = insertionBranch;  
  auto rightSubtree = insertionBranch->back;  
  std::set<unsigned int> leftGids;
  std::set<unsigned int> rightGids;
  // update distances between spid and the other nodes
  updateDistanceMatrixAux(gid, leftSubtree, 1, spidToGid, res, leftGids);
  updateDistanceMatrixAux(gid, rightSubtree, 1, spidToGid, res, rightGids);
  // update the distances between all the other nodes after adding spid
  for (auto gid1: leftGids) {
    for (auto gid2: rightGids) {
      res[gid1][gid2] += 1;
      res[gid2][gid1] += 1;
    }
  }
  return res;
}



void computeDiff(Node *currentInsertionNode,
    double lastDiff,
    Node *lastInsertionNode,
    const IntPairToDouble &deltas,
    std::vector<double> &diffVector,
    std::vector<double> &phi)

{
  double diff = 0.0;
  if (lastInsertionNode) {
    auto nodeB = lastInsertionNode;
    auto nodeC = currentInsertionNode->back;
    auto nodeA = StepwiseTree::getOtherNext(nodeB, nodeC);
    assert(nodeA->next && nodeB->next && nodeC->next);
    assert(nodeA->next == nodeB || nodeA->next == nodeC);
    std::pair<unsigned int, unsigned int> pairAB = {nodeA->index, nodeB->index};
    std::pair<unsigned int, unsigned int> pairAC = {nodeA->index, nodeC->index};
    // compute delta, the score different between currentInsertionNode
    // and lastInsertionNode
    double deltaAC = deltas.at(pairAC);
    double deltaAB = deltas.at(pairAB);
    double delta = 0.0;
    delta += 0.5 * (deltaAB - deltaAC);
    //std::cout << deltaAB << " " << deltaAC  << std::endl; 
    //std::cout << delta << std::endl;
    delta += 0.5 * (phi[nodeB->index] - phi[nodeC->index]);
    // add delta to the diff of the current node, and save it
    diff = lastDiff + delta;
    diffVector[currentInsertionNode->index] = diff;
    diffVector[currentInsertionNode->back->index] = diff;
    //std::cout << "diff = " << diff << std::endl ;
  }
  // recursion
  if (currentInsertionNode->next) {
    std::array<Node *, 2> children;
    children[0] = currentInsertionNode->next->back;
    children[1] =  currentInsertionNode->next->next->back;
    for (auto child: children) {
      computeDiff(child, 
          diff, 
          currentInsertionNode, 
          deltas,
          diffVector,
          phi);
    }
  }
  
}


Node *findBestInsertionBranch(unsigned int spid, // of the taxa to add
    StepwiseTree &tree,
    InducedStepwiseTrees &inducedTrees,
    std::vector<MatrixUint> &speciesMatrices,
    const std::vector<GeneCell *> &geneCells,
    UIntMatrix &spidToGid,
    std::vector<IntPairToDouble> &deltas,
    MatrixDouble &phis)
  
{
  Logger::info << "find best insertion " << spid << std::endl;
  auto bestScore = std::numeric_limits<double>::max();
  Node *bestBranch = nullptr;
  std::vector<double> perSuperBranchScore(tree.getAllNodes().size(), 0.0);

  for (unsigned int k = 0; k < geneCells.size(); ++k) {
    auto &inducedSpeciesTree = *inducedTrees[k];
    auto Nk = geneCells[k]->coverage.count();
    std::vector<double> diffVector(Nk * 4 - 6, 0.0);
    if (inducedSpeciesTree.getWrappedTree().getBranches().size() <= 2) {
      continue;
    }
    if (!geneCells[k]->coverage[spid]) {
      continue;
    }
    computeDiff(inducedSpeciesTree.getWrappedTree().getAnyLeaf()->back,
      0.0,
      nullptr,
      deltas[k],
      diffVector, 
      phis[k]);
    for (auto inducedBranch: inducedSpeciesTree.getWrappedTree().getBranches()) {
      for (auto superBranch: inducedSpeciesTree.getSuperBranches(inducedBranch)) {
        auto score = diffVector[inducedBranch->index];
        perSuperBranchScore[superBranch->index] += score;       
        perSuperBranchScore[superBranch->back->index] += score;       
      }
    }
  }
  for (auto branch: tree.getBranches()) {
    auto score = perSuperBranchScore[branch->index];
    ParallelContext::sumDouble(score);
    if (score < bestScore) {
      bestScore = score;
      bestBranch = branch;
    }
  }
  // now update the distance matrices
  for (unsigned int k = 0; k < geneCells.size(); ++k) {
    auto geneCell = geneCells[k];
    if (!geneCell->coverage[spid]) {
      continue;
    }
    auto &inducedSpeciesTree = *inducedTrees[k];
    auto &inducedSpeciesTreeMatrix = speciesMatrices[k];
    auto inducedBranch = inducedSpeciesTree.getInducedBranch(bestBranch);
    auto gid  = spidToGid[k][spid];
    inducedSpeciesTreeMatrix = getUpdatedDistanceMatrix(gid,
        inducedBranch,
        spidToGid[k],
        inducedSpeciesTreeMatrix);
  } 

  Logger::info << "Best score: " << bestScore << std::endl;
  return bestBranch;

}

void getSubClade(Node * node,
    std::set<Node *> &subclade)
{
  if (!node->next) {
    subclade.insert(node);
  } else {
    getSubClade(node->next->back, subclade);
    getSubClade(node->next->next->back, subclade);
  }
}

double StepwiseAsteroid::_computeDelta(std::set<Node *>&s1,
    std::set<Node *>s2,
    const DistanceMatrix &geneMatrix,
    const MatrixUint &speciesMatrix,
    const std::vector<unsigned int> &spidToGid)
{
  double delta = 0.0;
  for (auto node1: s1) {
    auto gid1 = spidToGid[node1->spid];
    for (auto node2: s2) {
      auto gid2 = spidToGid[node2->spid];
      delta += geneMatrix[gid1][gid2] * _pows[speciesMatrix[gid1][gid2]]; 
    }
  } 
  return delta;
}


double StepwiseAsteroid::_getPhi(unsigned int gid,
    Node *node,
    const DistanceMatrix &geneMatrix,
    const std::vector<unsigned int> &spidToGid,
    unsigned int depth)
{
  if (!node->next) {
    auto gid2 = spidToGid[node->spid];
    return geneMatrix[gid][gid2] * _pows[depth];
  } else {
    double res = 0.0;
    res += _getPhi(gid, 
        node->next->back, 
        geneMatrix, 
        spidToGid, 
        depth + 1);
    res += _getPhi(gid, 
        node->next->next->back, 
        geneMatrix, 
        spidToGid, 
        depth + 1);
    return res;
  }
}

/*
void fillDeltaForAddedTaxon(Node *node,
    unsigned int taxonGid,
    const DistanceMatrix &geneDistances,
    DistanceMatrix &deltas)
{
  if (!node->next) {
  }
}
*/

void StepwiseAsteroid::_updatePhis(unsigned int spid)
{
  for (unsigned int k = 0; k < _K; ++k) {
    auto gid = _spidToGid[k][spid];
    for (auto node: _inducedTrees[k]->getWrappedTree().getAllNodes()) {
      _phis[k][node->index] = _getPhi(gid,
          node, 
          _geneCells[k]->distanceMatrix,
          _spidToGid[k],
          0);
    }
  }
}
      
void StepwiseAsteroid::_updateDeltas()
{
  for (unsigned int k = 0; k < _K; ++k) {
    for (auto node: _inducedTrees[k]->getWrappedTree().getAllNodes()) {
      if (!node->next) {
        continue;
      }
      // blabla
      std::set<Node *> cladeB;
      std::set<Node *> cladeC;
      auto nodeB = node->next->back;
      auto nodeC = node->next->next->back;
      getSubClade(nodeB, cladeB);
      getSubClade(nodeC, cladeC);
      auto d = _computeDelta(cladeB,
          cladeC,
          _geneCells[k]->distanceMatrix,
          _speciesMatrices[k],
          _spidToGid[k]);
      std::pair<unsigned int, unsigned int> p1, p2;
      p1 = {nodeB->back->index, nodeC->back->index};
      p2 = {nodeC->back->index, nodeB->back->index};
      _slowDeltas[k][p1] = d;
      _slowDeltas[k][p2] = d;
    }

  }
}

void StepwiseAsteroid::insertLabel(const std::string &label)
{
  auto spid = _speciesToSpid.at(label);
  Logger::info << "Adding " << label  << " spid=" << spid << std::endl;
  _insertedSpecies.set(spid);
  if (_insertedSpecies.count() == 1) {
    _tree.addLeaf(spid, nullptr);
    return;
  } else if (_insertedSpecies.count() == 2) {
    _tree.addLeaf(spid, _tree.getAnyBranch());
    return;
  }
  _updateDeltas();
  _updatePhis(spid);
  auto bestBranch = findBestInsertionBranch(spid,
      _tree,
      _inducedTrees,
      _speciesMatrices,
      _geneCells,
      _spidToGid,
      _slowDeltas,
      _phis);
  _tree.addLeaf(spid, bestBranch);
}

void StepwiseAsteroid::exportTree(const std::string &outputPath)
{
  ParallelOfstream os(outputPath);
  os << _tree.getNewickString(_spidToSpecies) << std::endl;

}

