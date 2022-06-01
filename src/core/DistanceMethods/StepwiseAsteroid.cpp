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
  _deltas(_K)

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
    _deltas[k] = MatrixDouble(internalNodes, 
        std::vector<double>(internalNodes, 0.0));
    auto &spidToGid = _spidToGid[k];
    auto &gidToSpid = _geneCells[k]->gidToSpid;
    spidToGid.resize(_N);
    for (unsigned int gid = 0; gid < gidToSpid.size(); ++gid) {
      spidToGid[gidToSpid[gid]] = gid;
    }
  }
}
  

static void computeDiffNew(Node *currentInsertionNode,
    double lastDiff,
    Node *lastInsertionNode,
    unsigned int newNodeIndex,
    const MatrixDouble &deltas,
    std::vector<double> &diffVector)
{
  double diff = 0.0;
  if (lastInsertionNode) {
    auto nodeC = lastInsertionNode;
    auto nodeB = currentInsertionNode->back;
    auto nodeA = StepwiseTree::getOtherNext(nodeB, nodeC);
    assert(nodeA->next && nodeB->next && nodeC->next);
    assert(nodeA->next == nodeB || nodeA->next == nodeC);
  
    nodeA = nodeA->back;
    nodeB = nodeB->back;
    nodeC = nodeC->back;
    // compute delta, the score different between currentInsertionNode
    // and lastInsertionNode
    double deltaAB = deltas[nodeA->index][nodeB->index];
    double deltaAC = deltas[nodeA->index][nodeC->index];
    auto t = newNodeIndex;
   // currentInsertionNode->index;
    double deltaBt = deltas[t][nodeB->index];
    double deltaCt = deltas[t][nodeC->index];
    
    assert(deltaBt * deltaCt != 0.0);
    assert (deltaAB * deltaAC != 0.0);
    double delta = 0.0;
    delta += 0.5 * (deltaAC - deltaAB);
    delta += 0.5 * (deltaBt - deltaCt);
    // add delta to the diff of the current node, and save it
    diff = lastDiff + delta;
    //std::cout << "  diff " << diff << std::endl;
    diffVector[currentInsertionNode->index] = diff;
    diffVector[currentInsertionNode->back->index] = diff;
  }
  // recursion
  if (currentInsertionNode->next) {
    std::array<Node *, 2> children;
    children[0] = currentInsertionNode->next->back;
    children[1] =  currentInsertionNode->next->next->back;
    for (auto child: children) {
      computeDiffNew(child, 
          diff, 
          currentInsertionNode, 
          newNodeIndex,
          deltas,
          diffVector);
    }
  }
  
}


static Node *findBestInsertionBranch(unsigned int spid, // of the taxa to add
    StepwiseTree &tree,
    InducedStepwiseTrees &inducedTrees,
    const std::vector<GeneCell *> &geneCells,
    const std::vector<MatrixDouble> &newDeltas)
  
{
  Logger::info << "find best insertion " << spid << std::endl;
  auto bestScore = std::numeric_limits<double>::max();
  Node *bestBranch = nullptr;
  std::vector<double> perSuperBranchScore(tree.getAllNodes().size(), 0.0);
  for (unsigned int k = 0; k < geneCells.size(); ++k) {
    auto &inducedSpeciesTree = *inducedTrees[k];
    auto Nk = geneCells[k]->coverage.count();
    std::vector<double> diffVector(Nk * 4 - 6, 0.0);
    std::vector<double> diffVectorNew(Nk * 4 - 6, 0.0);
    if (inducedSpeciesTree.getWrappedTree().getBranches().size() <= 2) {
      continue;
    }
    if (!geneCells[k]->coverage[spid]) {
      continue;
    }
    auto taxonIndex = inducedTrees[k]->getWrappedTree().getNextNewNodeIndex();
    computeDiffNew(inducedSpeciesTree.getWrappedTree().getAnyLeaf()->back,
      0.0,
      nullptr,
      taxonIndex,
      newDeltas[k],
      diffVectorNew);
    for (auto inducedBranch: inducedSpeciesTree.getWrappedTree().getBranches()) {
      
      for (auto superBranch: inducedSpeciesTree.getSuperBranches(inducedBranch)) {
        auto score = diffVectorNew[inducedBranch->index];
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
  //Logger::info << "Best score: " << bestScore << std::endl;
  return bestBranch;

}



double StepwiseAsteroid::_updateDeltasNewTaxonAux(
    unsigned int taxonGid,
    unsigned int taxonIndex,
    unsigned int k,
    Node *node)
{
  double d = 0.0;
  if (!node->next) { // leaf case
    auto gid = _spidToGid[k][node->spid];
    d = _geneCells[k]->distanceMatrix[taxonGid][gid]; 
    _deltas[k][taxonIndex][node->index] = d; 
    _deltas[k][node->index][taxonIndex] = d;
    return d;
  } else {
    auto left = node->next->back;
    auto right = node->next->next->back;
    d += _updateDeltasNewTaxonAux(taxonGid, taxonIndex, k, left);
    d += _updateDeltasNewTaxonAux(taxonGid, taxonIndex, k, right);
    d *= 0.5;
  }
  return d;
}

void StepwiseAsteroid::_updateDeltasNewTaxon(unsigned int taxonSpid)
{
  for (unsigned int k = 0; k < _K; ++k) {
    if (!_geneCells[k]->coverage[taxonSpid]) {
      continue;
    }
    auto anyNode = _inducedTrees[k]->getAnyNode();
    if (!anyNode) {
      // the induced tree is empty, nothing to do
      continue;
    }
    // get the index of the node to be added. The node does not exist yet
    // but we know which index will be assigned to it
    auto taxonIndex = _inducedTrees[k]->getWrappedTree().getNextNewNodeIndex();
    auto taxonGid = _spidToGid[k][taxonSpid];
    std::vector<Node *> postOrderNodes;
    _inducedTrees[k]->getWrappedTree().fillPostOrder(postOrderNodes);
    for (auto node: postOrderNodes) {
      if (!node->next) { // leaf case
        auto nodeGid = _spidToGid[k][node->spid];
        auto d = _geneCells[k]->distanceMatrix[taxonGid][nodeGid]; 
        _deltas[k][taxonIndex][node->index] = d; 
        _deltas[k][node->index][taxonIndex] = d;
        assert(d != 0.0);
      } else {
        auto left = node->next->back;
        auto right = node->next->next->back;
        auto d = 0.5 *
          (_deltas[k][left->index][taxonIndex] + _deltas[k][right->index][taxonIndex]);
        _deltas[k][taxonIndex][node->index] = d; 
        _deltas[k][node->index][taxonIndex] = d;
        assert(d != 0.0);
      }
    }
  }
}

/*
  insertedNode   Y1              Yn
        |         |               |
        *         *               *
  *----* *-------* *---- .... ---* *---------* 
  ^      ^         ^               ^         ^
  A  chain[0]    chain[1]     chain.back()  B

*/
void StepwiseAsteroid::_updateDeltasChain(unsigned int k,
    Node *insertedNode,
    const std::vector<Node *> &chain)
{
  //auto A = chain[0]->back;
  const auto A = StepwiseTree::getOtherNext(insertedNode->back, chain[0])->back; 
  const auto B = chain.back()->back;
  for (unsigned int i = 0; i < chain.size(); ++i) {
    auto C = chain[i];
    auto Cminus1 = i ? chain[i-1] : A;
    auto Yback = StepwiseTree::getOtherNext(C, Cminus1->back);
    assert(Yback);
    auto Y = Yback->back;
    assert(Y);
    auto d = 0.5 * 
    (_deltas[k][B->index][Cminus1->index]
     + _deltas[k][B->index][Y->index]);
    _deltas[k][B->index][C->index] = d;
    _deltas[k][C->index][B->index] = d;
    assert (d != 0.0);
  }
}

void StepwiseAsteroid::_updateDeltasAux(unsigned int k,
    Node *insertedNode,
    Node *Bback,
    std::vector<Node *> &chain)
{
  chain.push_back(Bback);
  _updateDeltasChain(k, insertedNode, chain);
  if (Bback->back->next) {
    auto Bback1 = Bback->back->next;//B->next->back;
    auto Bback2 = Bback->back->next->next;
    //auto B2 = B->next->next->back;
    _updateDeltasAux(k, insertedNode, Bback1, chain);
    _updateDeltasAux(k, insertedNode, Bback2, chain);
  }
  chain.pop_back();
}

void StepwiseAsteroid::_updateDeltas()
{
  for (unsigned int k = 0; k < _K; ++k) { // for each induced tree
    // inserted node is the node that was inserted in the induced
    // tree k at the last step 
    auto insertedNode = _tree.getLastAddedInducedNodes()[k];
    if (!insertedNode) {
      // this induced tree does not cover the taxon inserted
      // at the last step
      continue;
    }
    // We now want to iterate over all nodes A that have the
    // inserted node as descendant, and for each such A, over
    // each node B whose subtree does not intersect the subtree 
    // of A. We also need to keep track of the chain chainAB of nodes
    // between the A and B nodes. We can then call _updateDeltasChain
    // on each such A and each chainAB
    //
    if (insertedNode->back && insertedNode->back->next) {
      std::vector<Node *> chain;
      auto Bback1 = insertedNode->back->next;
      auto Bback2 = insertedNode->back->next->next;
      _updateDeltasAux(k, insertedNode, Bback1, chain);
      _updateDeltasAux(k, insertedNode, Bback2, chain);
      // we need to update one last avaerage distance: the distance 
      // between the new taxon and its back node that was just inserted
      // and thus not updated with _updateDeltasNewTaxon
      auto d = _deltas[k][insertedNode->index][Bback1->back->index];
      d += _deltas[k][insertedNode->index][Bback2->back->index];
      d *= 0.5;
      _deltas[k][insertedNode->index][insertedNode->back->index] = d;
      _deltas[k][insertedNode->back->index][insertedNode->index] = d;

    }
  }
   
}


void StepwiseAsteroid::insertLabel(const std::string &label)
{
  auto spid = _speciesToSpid.at(label);
  Logger::info << "Adding " << label  << " spid=" << spid << std::endl;
  _insertedSpecies.set(spid);
  if (_insertedSpecies.count() == 1) {
    // this is the first taxon to add
    _tree.addLeaf(spid, nullptr);
    return;
  } else if (_insertedSpecies.count() == 2) {
    // this is the second taxon to add
    _updateDeltasNewTaxon(spid); 
    _tree.addLeaf(spid, _tree.getAnyBranch());
    return;
  }
  _updateDeltasNewTaxon(spid); 
  auto bestBranch = findBestInsertionBranch(spid,
      _tree,
      _inducedTrees,
      _geneCells,
      _deltas);
  _tree.addLeaf(spid, bestBranch);
  _updateDeltas();
}

void StepwiseAsteroid::exportTree(const std::string &outputPath)
{
  ParallelOfstream os(outputPath);
  os << _tree.getNewickString(_spidToSpecies) << std::endl;
}

