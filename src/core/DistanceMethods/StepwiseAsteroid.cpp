#include "StepwiseAsteroid.hpp"

#include <limits>
#include <IO/ParallelOfstream.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <set>

StepwiseAsteroid::StepwiseAsteroid(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells):
  _N(speciesToSpid.size()),
  _K(geneCells.size()),
  _speciesToSpid(speciesToSpid),
  _spidToSpecies(_N),
  _geneCells(geneCells),
  _insertedSpecies(_N, false),
  _speciesMatrices(_K)

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
    _speciesMatrices[k] = DistanceMatrix(Nk, std::vector<double>(Nk, 0.0));

  }
}
  


double getAsteroidScore(const DistanceMatrix &speciesMatrix,
    const std::vector<GeneCell*> &geneCells,
    const BitVector &addedSpids) // the spids that belong to the tree being built
{
  double score = 0.0;
  for (auto geneCell: geneCells) {
    if ((addedSpids & geneCell->coverage).count() < 4) {
      continue;
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
        double m = speciesMatrix[spid1][spid2];
        score += d * pow(0.5, m);
      }
    }
  }
  //ParallelContext::sumDouble(score);
  return score;
}

/**
 *  Compute the internode distance bitween spid and the leaves under node
 *  Also fill subspids with the spid of the leaves under node
 */
void updateDistanceMatrixAux(unsigned int spid,
    Node *node,
    double  distance,
    DistanceMatrix &distanceMatrix,
    std::set<unsigned int> &subspids)
{
  if (!node->next) {
    auto spid2 = node->spid;
    distanceMatrix[spid][spid2] = distance;
    distanceMatrix[spid2][spid] = distance;
    subspids.insert(spid2);
    return;
  }
  auto left = node->next->back;
  auto right = node->next->next->back;
  updateDistanceMatrixAux(spid, left, distance + 1.0, distanceMatrix, subspids);
  updateDistanceMatrixAux(spid, right, distance + 1.0, distanceMatrix, subspids);
}
    

void printMatrix(const DistanceMatrix &m) {
  for (const auto &v: m) {
    for (const auto val: v) {
      Logger::info << val << " ";
    }
    Logger::info << std::endl;
  }
}

DistanceMatrix getUpdatedDistanceMatrix(unsigned int spid,
    Node *insertionBranch,
    const DistanceMatrix &speciesMatrix
    )
{
  DistanceMatrix res(speciesMatrix);
  if (!insertionBranch || !insertionBranch->back) {
    // edge case:
    // this is the second inserted taxa
    // we don't need to update the matrix
    return res;
  }
  auto leftSubtree = insertionBranch;  
  auto rightSubtree = insertionBranch->back;  
  std::set<unsigned int> leftSpids;
  std::set<unsigned int> rightSpids;
  // update distances between spid and the other nodes
  updateDistanceMatrixAux(spid, leftSubtree, 1.0, res, leftSpids);
  updateDistanceMatrixAux(spid, rightSubtree, 1.0, res, rightSpids);
  // update the distances between all the other nodes after adding spid
  for (auto spid1: leftSpids) {
    for (auto spid2: rightSpids) {
      res[spid1][spid2] += 1.0;
      res[spid2][spid1] += 1.0;
    }
  }
  return res;
}




Node *findBestInsertionBranch(unsigned int spid, // of the taxa to add
    StepwiseTree &tree,
    InducedStepwiseTrees &inducedTrees,
    std::vector<DistanceMatrix> &speciesMatrices,
    const std::vector<GeneCell *> &geneCells,
    const BitVector &addedSpids) // the spids that belong to the tree being built
{
  Logger::info << "find best insertion " << spid << std::endl;
  auto bestScore = std::numeric_limits<double>::max();
  Node *bestBranch = nullptr;
  std::vector<double> perSuperBranchScore(tree.getAllNodes().size(), 0.0);
  
  for (unsigned int k = 0; k < geneCells.size(); ++k) {
    auto geneCell = geneCells[k];
    if (!geneCell->coverage[spid]) {
      // spid is not covered by this gene tree and thus
      // does not contribute to the score diff
      continue;
    }
    auto &inducedSpeciesTree = *inducedTrees[k];
    auto &inducedSpeciesTreeMatrix = speciesMatrices[k];
    for (auto inducedBranch: inducedSpeciesTree.getWrappedTree().getBranches()) {
      auto newSpeciesMatrix = getUpdatedDistanceMatrix(spid,
          inducedBranch,
          inducedSpeciesTreeMatrix);
      auto score = getAsteroidScore(newSpeciesMatrix,
          geneCells,
          addedSpids);
      for (auto superBranch: inducedSpeciesTree.getSuperBranches(inducedBranch)) {
        perSuperBranchScore[superBranch->index] += score;       
        perSuperBranchScore[superBranch->back->index] += score;       
      }
    }
  }
  for (auto branch: tree.getBranches()) {
    auto score = perSuperBranchScore[branch->index];
    ParallelContext::sumDouble(score);
    //std::cout << score << std::endl;
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
    inducedSpeciesTreeMatrix = getUpdatedDistanceMatrix(spid,
        inducedBranch,
        inducedSpeciesTreeMatrix);
  
  }  
  Logger::info << "Best score: " << bestScore << std::endl;
  // TODO:
  // - for each induced tree, for each of its insertion branches: 
  //   * compute its score diff
  //   * update the relevant values in perSuperBranchDiff 


  /*
  for (auto branch: tree.getBranches()) {
    auto newSpeciesMatrix = getUpdatedDistanceMatrix(spid, 
        branch,
        speciesMatrix);
    auto score = getAsteroidScore(newSpeciesMatrix,
        geneCells,
        addedSpids);
    if (score < bestScore) {
      bestScore = score;
      bestBranch = branch;
    }
  }
  Logger::info << "Best score " << bestScore << std::endl;
  if (bestBranch) {
    speciesMatrix = getUpdatedDistanceMatrix(spid, bestBranch, speciesMatrix);
  }
  */
  return bestBranch;

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
  auto bestBranch = findBestInsertionBranch(spid,
      _tree,
      _inducedTrees,
      _speciesMatrices,
      _geneCells,
      _insertedSpecies);
  _tree.addLeaf(spid, bestBranch);
}

void StepwiseAsteroid::exportTree(const std::string &outputPath)
{
  ParallelOfstream os(outputPath);
  os << _tree.getNewickString(_spidToSpecies) << std::endl;

}

