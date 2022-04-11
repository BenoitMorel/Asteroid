#include "StepwiseAsteroid.hpp"

#include <limits>
#include <IO/ParallelOfstream.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <set>

StepwiseAsteroid::StepwiseAsteroid(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells,
    const std::array<std::string, 3> initialLabels):
  _N(speciesToSpid.size()),
  _K(geneCells.size()),
  _speciesToSpid(speciesToSpid),
  _spidToSpecies(_N),
  _geneCells(geneCells),
  _insertedSpecies(_N, false),
  _speciesMatrix(_N, std::vector<double>(_N, 0.0))

{
  for (const auto &label1: initialLabels) {
    _insertedSpecies.set(_speciesToSpid.at(label1));
    auto spid1 = speciesToSpid.at(label1);
    _tree.addLeaf(spid1, _tree.getAnyBranch());
    for (const auto &label2: initialLabels) {
      if (label2 != label1) {
        auto spid2 = speciesToSpid.at(label2);
        _speciesMatrix[spid1][spid2] = 1.0;
      }
    }
  }
  for (auto pair: speciesToSpid) {
    auto label = pair.first;
    auto spid = pair.second;
    _spidToSpecies[spid] = label;
  }
}
  


double getAsteroidScore(const DistanceMatrix &speciesMatrix,
    const std::vector<GeneCell*> &geneCells,
    const BitVector &addedSpids) // the spids that belong to the tree being built
{
  double score = 0.0;
  for (auto geneCell: geneCells) {
    auto inter = addedSpids & geneCell->coverage;
    if (inter.count() < 4) {
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
  ParallelContext::sumDouble(score);
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
    

DistanceMatrix getUpdatedDistanceMatrix(unsigned int spid,
    Node *insertionBranch,
    const DistanceMatrix &speciesMatrix
    )
{
  DistanceMatrix res(speciesMatrix);
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



void printMatrix(const DistanceMatrix &m) {
  for (const auto &v: m) {
    for (const auto val: v) {
      Logger::info << val << " ";
    }
    Logger::info << std::endl;
  }
}

Node *findBestInsertionBranch(unsigned int spid, // of the taxa to add
    StepwiseTree &tree,
    DistanceMatrix &speciesMatrix,
    const std::vector<GeneCell *> &geneCells,
    const BitVector &addedSpids) // the spids that belong to the tree being built
{
  auto bestScore = std::numeric_limits<double>::max();
  Node *bestBranch = nullptr;
  for (auto branch: tree.getBranches()) {
    auto newSpeciesMatrix = getUpdatedDistanceMatrix(spid, 
        branch,
        speciesMatrix);
    //Logger::info << "New distance matrix: " << std::endl;
    //printMatrix(newSpeciesMatrix);
    auto score = getAsteroidScore(newSpeciesMatrix,
        geneCells,
        addedSpids);
    if (score < bestScore) {
      bestScore = score;
      bestBranch = branch;
    }
  }
  Logger::info << "Best score " << bestScore << std::endl;
  speciesMatrix = getUpdatedDistanceMatrix(spid, bestBranch, speciesMatrix);
  assert(bestBranch);
  return bestBranch;

}

void StepwiseAsteroid::insertLabel(const std::string &label)
{
  auto spid = _speciesToSpid.at(label);
  Logger::info << "Adding " << label  << " spid=" << spid << std::endl;
  _insertedSpecies.set(spid);
  auto bestBranch = findBestInsertionBranch(spid,
      _tree,
      _speciesMatrix,
      _geneCells,
      _insertedSpecies);
  _tree.addLeaf(spid, bestBranch);
}

void StepwiseAsteroid::exportTree(const std::string &outputPath)
{
  ParallelOfstream os(outputPath);
  os << _tree.getNewickString(_spidToSpecies) << std::endl;

}

