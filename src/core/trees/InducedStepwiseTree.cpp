#include "InducedStepwiseTree.hpp"
#include <iostream>
  
InducedStepwiseTree::InducedStepwiseTree(const StepwiseTree &superTree,
      const BitVector &coverage): _supertree(superTree),
  _coverage(coverage),
  _superToInduced(4 * coverage.size() - 6),
  _inducedToSuper(4 * coverage.count() - 6),
  _inducedToSuperGhost(4 * coverage.count() - 6),
  _taxaCount(0)
{
}

void InducedStepwiseTree::mapNodes(Node *superNode,
    Node *inducedNode)
{
  std::cout << "map node  " << superNode->index << std::endl;
  assert(_inducedToSuper[inducedNode->index] == nullptr);
  //assert(_superToInduced[superNode->index]  == nullptr);
  _inducedToSuper[inducedNode->index] = superNode;
  _superToInduced[superNode->index] = inducedNode;
}
  
void InducedStepwiseTree::mapGhostNode(Node *superNode,
      Node *inducedBranch)
{
  std::cout << "map ghost  " << superNode->index << std::endl;
  //assert(nullptr ==  _superToInduced[superNode->index]);
  _inducedToSuperGhost[inducedBranch->index].insert(superNode);
  _inducedToSuperGhost[inducedBranch->back->index].insert(superNode);
  _superToInduced[superNode->index] = inducedBranch;
}


void InducedStepwiseTree::unmapGhostNode(Node *superNode)
{
  std::cout << "unmap node  " << superNode->index << std::endl;;
  auto inducedBranch = _superToInduced[superNode->index];
  //assert(inducedBranch);
  if (!inducedBranch) {
    // already unmapped
    return;
  }
  _superToInduced[superNode->index] = nullptr;
  _inducedToSuperGhost[inducedBranch->index].erase(superNode);
  _inducedToSuperGhost[inducedBranch->back->index].erase(superNode);
}

static bool areSameInnerNode(Node *n1, Node *n2)
{
  assert(n1 != n2);
  return n1->next == n2 || n1->next->next == n2;
}

bool findSuperLCAs(Node *node, 
        Node *superBound2,
        Node *superBound3,
        bool &superBound2IsInSubtree,
        bool &superBound3IsInSubtree,
        std::array<Node *, 3> &lcas)
{
  std::cout <<"  treating node " << node->index << std::endl;
  if (node == superBound2) {
    lcas[1] = node;
    std::cout << "found 2 " << std::endl;
    superBound2IsInSubtree = true;
    return true;
  } else if (node == superBound3) {
    lcas[2] = node;
    superBound3IsInSubtree = true;
    std::cout << "found 3 " << std::endl;
    return true;
  }
  if (!node->next) {
    return false;
  }

  auto left = node->next->back;
  auto right = node->next->next->back;
  bool superBound2IsInSubtreeLeft = false;
  bool superBound2IsInSubtreeRight = false;
  bool superBound3IsInSubtreeLeft = false;
  bool superBound3IsInSubtreeRight = false;
  //std::cout << "go down left" << std::endl;
  auto foundLeft = findSuperLCAs(left, 
      superBound2, 
      superBound3, 
      superBound2IsInSubtreeLeft,
      superBound3IsInSubtreeLeft,
      lcas);
  
  //std::cout << "go up left" << std::endl;
  if (lcas[0]) { // we have found our lcas, exit
    std::cout << "exit recursion" << std::endl;
    return true;
  } 
  //std::cout << "go down right" << std::endl;
  auto foundRight = findSuperLCAs(right, 
      superBound2, 
      superBound3, 
      superBound2IsInSubtreeRight,
      superBound3IsInSubtreeRight,
      lcas);
  superBound2IsInSubtree = superBound2IsInSubtreeLeft || 
    superBound2IsInSubtreeRight;
  superBound3IsInSubtree = superBound3IsInSubtreeLeft || 
    superBound3IsInSubtreeRight;
  //std::cout << "go up right" << std::endl;
  if (lcas[0]) { // we have found our lcas, exit
    std::cout << "exit recursion" << std::endl;
    return true;
  }
  if (foundLeft && foundRight) {
    std::cout << "we are at the LCA!" << std::endl;
    assert(superBound2IsInSubtree && superBound3IsInSubtree);
    // we are at the LCA!
    // store the final lca value
    // now that lcas[0] != nullptr, all recursive
    // calls will be skipped
    lcas[0] = node;
    lcas[1] = lcas[1]->back;
    lcas[2] = lcas[2]->back;
    assert(areSameInnerNode(lcas[0], lcas[1]));
    assert(areSameInnerNode(lcas[0], lcas[2]));
    assert(areSameInnerNode(lcas[1], lcas[2]));
    return true;
  }
  // if one of the subtrees contains either superBound2
  // or superBound3, we update the corresponding cell
  // to node: lcas[1] is the subtree countaining superBound2
  // and lcas[2] is the subtree countaining superBound3 (if
  // they were found)
  if (superBound2IsInSubtree) {
    assert(!superBound3IsInSubtree);
    std::cout << "up superbound2" << std::endl;
    lcas[1] = node;
  }
  if (superBound3IsInSubtree) {
    assert(!superBound2IsInSubtree);
    std::cout << "up superbound3" << std::endl;
    lcas[2] = node;
  }
  return foundRight || foundLeft;
  
}
 

void InducedStepwiseTree::addLeaf(Node *superLeaf,
      unsigned int spid)
{ 
  
  
  std::cout << "current supertree:" << std::endl;
  std::vector<std::string> plop; 
  std::cout << _supertree.getNewickString(plop) << std::endl;
  // this function has many special cases to consider
  // - is the new spid covered by this induced gene tree?
  // - how many taxa do we already have in the current induced tree? (0, 1, or more)

  if (_taxaCount == 0) {
    // first edge case, the induced tree is currently empty
    if (_coverage[spid]) {
      // if the taxa is covered, add it to the induced tree
      auto inducedLeaf = _inducedTree.addLeaf(spid, nullptr);
      _taxaCount++;
      // and map it to its corresponding leaf in the supertree
      mapNodes(superLeaf, inducedLeaf);
    } 
    // if the taxa is not covered, we do not have anything to do
    return;
  } 
  if (_taxaCount == 1) {
    // second edge case: we only have one leaf
    if (_coverage[spid]) {
      // if the taxa is covered, add it to the induced tree
      auto inducedLeaf = _inducedTree.addLeaf(spid, _inducedTree.getAnyBranch());
      _taxaCount++;
      // and map it to its corresponding leaf in the supertree
      mapNodes(superLeaf, inducedLeaf);
      // get the second leaf that is already in the induced tree
      auto otherInducedLeaf = inducedLeaf->back;
      auto otherSuperLeaf = _inducedToSuper[otherInducedLeaf->index];
      // now we have one single induced branch: we map to this branch 
      // all the nodes of the supertree that are not the two leaves
      // covered by the current induced tree
      auto singleInducedBranch = inducedLeaf;
      for (auto superNode: _supertree.getAllNodes()) {
        if (superNode != superLeaf && superNode != otherSuperLeaf) {
          mapGhostNode(superNode, singleInducedBranch);
        }
      }
    }
    // if the taxa is not covered, we do not have anything to do
    return;
  }
  
  int unmappedCount = 0;
  for (auto node: _supertree.getAllNodes()) {
    if (!_superToInduced[node->index]) {
      unmappedCount++;
      std::cout << "unmapped: " << node->index << std::endl;
    }
  }
  for (auto node: _inducedTree.getAllNodes()) {
    assert(_inducedToSuper[node->index]);
  }
  assert(unmappedCount == 4);
  
  // the three freshly added internal nodes in the supertree:
  auto super1 = superLeaf->back;
  auto super2 = super1->next;
  auto super3 = super2->next;
  if (!_coverage[spid]) {
    // the induced tree does not cover this taxa
    // so we do not add any node to the induced tree
    // but we still have to map the new super nodes
    // as "ghost nodes" to the insertion branch in the induced tree 
    auto inducedBranch = _superToInduced[super2->back->index];
    mapGhostNode(super1, inducedBranch); 
    mapGhostNode(super2, inducedBranch); 
    mapGhostNode(super3, inducedBranch); 
    mapGhostNode(superLeaf, inducedBranch);
  } else {
    // the induced tree covers this taxa, so we add it
    // to the induced tree
    
    // insertionBranch: the branch at which we insert the new
    // taxa in the induced tree
    auto insertionBranch = _superToInduced[super2->back->index];
    auto insertionBranchGhosts = _inducedToSuperGhost[insertionBranch->index];
    // the super nodes mapped to the added leaf 
    // and to the two extremities of the insertion branch. 
    // We need to get them before inserting the induced leaf
    auto superBound1 = superLeaf;
    auto superBound2 = _inducedToSuper[insertionBranch->index];
    auto superBound3 = _inducedToSuper[insertionBranch->back->index];
    // we unmap all ghost nodes mapped to the insertion branch
    // because we will reassign all of them to one of
    // the 3 branches starting from the LCAs
    for (auto ghost: insertionBranchGhosts) {
      unmapGhostNode(ghost);
    }
    // insert the leaf into the induced tree    
    auto inducedLeaf = _inducedTree.addLeaf(spid, 
        insertionBranch);
    std::cout << "super leaf: " << superBound1->index << std::endl;
    // the three freshly added internal nodes
    NodeTriplet inducedTriplet;
    inducedTriplet[0] = inducedLeaf->back;
    inducedTriplet[1] = inducedTriplet[0]->next;
    inducedTriplet[2] = inducedTriplet[1]->next;
    

    // now we need to find the internal node that separates superBound1/2/3
    std::array<Node *, 3> lcas;
    lcas[0] = lcas[1] = lcas[2] = nullptr;
    bool dummy1 = false;
    bool dummy2 = false;
    std::cout << "looking for superbound2=" << superBound2->index << " and superbound3=" << superBound3->index << "from " << superBound1->index <<  std::endl;
    findSuperLCAs(superBound1->back, 
        superBound2,
        superBound3,
        dummy1,
        dummy2,
        lcas);
    assert(lcas[0]);
    assert(lcas[1]);
    assert(lcas[2]);
    std::cout << "map leaf to superleaf " << superLeaf->index << std::endl;
    mapNodes(superLeaf, inducedLeaf);
    assert(_superToInduced[superBound2->index] == inducedTriplet[1]->back);
    for (unsigned int i = 0; i < 3; ++i) {
      assert(lcas[i]);
      std::cout << "map LCA: " ;
      mapNodes(lcas[i], inducedTriplet[i]);
    }
    std::cout << "Map ghosts until bound 1" << std::endl;
    mapGhostUntilBound(lcas[0]->back, superBound1, inducedTriplet[0]);
    std::cout << "Map ghosts until bound 2" << std::endl;
    mapGhostUntilBound(lcas[1]->back, superBound2, inducedTriplet[1]);
    std::cout << "Map ghosts until bound 3" << std::endl;
    mapGhostUntilBound(lcas[2]->back, superBound3, inducedTriplet[2]);
    /*
    for (unsigned int i = 0; i < 3; ++i) {
      for (auto superGhost: tripartition[i]) {
        mapGhostNode(superGhost, inducedTriplet[i]);
      }
    }
    insertionBranchGhosts.insert(super2);
    insertionBranchGhosts.insert(super3);
    for (auto ghost: insertionBranchGhosts) {
      if (!_superToInduced[ghost->index]) {
        mapGhostNode(ghost, inducedTriplet[0]);
      }
    }
    */
  }
}

void InducedStepwiseTree::mapGhostUntilBound(Node *superGhost,
      Node *superBound,
      Node *inducedToBeAttached,
      bool firstCall)
{
  if (!firstCall) { 
    mapGhostNode(superGhost->back, inducedToBeAttached);
  }
  if (superGhost == superBound) {
    return;
  }
  mapGhostNode(superGhost, inducedToBeAttached);
  if (superGhost->next) {
    mapGhostUntilBound(superGhost->next->back,
        superBound,
        inducedToBeAttached,
        false);
    mapGhostUntilBound(superGhost->next->next->back,
        superBound,
        inducedToBeAttached,
        false);
  }

}

