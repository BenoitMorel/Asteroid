#pragma once

#include <unordered_set>
#include <vector>
#include <maths/bitvector.hpp>
#include <trees/StepwiseTree.hpp>

using BranchSet = std::unordered_set<Node *>;
using InducedToSuperNodes = std::vector<Node *>;
using SuperToInducedNodes = std::vector<Node *>;
using InducedToSuperGhost = std::vector<BranchSet>;
using NodeTriplet = std::array<Node *, 3>;
class InducedStepwiseTree
{
public:
  InducedStepwiseTree(const StepwiseTree &superTree,
      const BitVector &coverage);
  virtual ~InducedStepwiseTree() {}

  void addLeaf(Node *superLeaf,
      unsigned int spid);

  const StepwiseTree &getWrappedTree() const {return _inducedTree;}

  BranchSet getSuperBranches(Node *inducedNode) const;
  Node *getInducedBranch(Node *superNode) const {return _superToInduced[superNode->index];}

private:
  // the induced tree beeing built
  StepwiseTree _inducedTree;
  // the supertree from which we build the induced tree
  const StepwiseTree &_supertree;
  // _coverage[spid] tells if spid should be covered 
  // by this induced tree
  BitVector _coverage;
  
  // mapping from a super node to the corresponding
  // induced node. If the super node is a "ghost" node
  // then it points to one extremity of the closest 
  // induced branch
  SuperToInducedNodes _superToInduced;
  // mapping 
  InducedToSuperNodes _inducedToSuper;
  InducedToSuperGhost _inducedToSuperGhost;

  // number of taxa that have already been added
  // (sometimes include the taxa that we are trying to add)
  unsigned int _taxaCount;


  // map superNode to inducedNode. The superNode must should
  // not be a ghost node.
  void mapNodes(Node *superNode,
    Node *inducedNode);

  // map superNode to both extremities of inducedBranch
  // as ghost nodes (nodes that do not belong to the
  // induced tree but that are assigned to the closest branch)
  void mapGhostNode(Node *superNode,
      Node *inducedBranch);

  // undo mapGhostNode
  void unmapGhostNode(Node *superNode);

  void mapGhostUntilBound(Node *superGhost,
      Node *superBound,
      Node *inducedToBeAttached,
      bool firstCall = true);

};


