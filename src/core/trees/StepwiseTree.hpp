#pragma once

#include <string>
#include <vector>

class InducedStepwiseTree;

struct Node {
  Node(): index(static_cast<unsigned int>(-1)), 
  back(nullptr), 
  next(nullptr), 
  spid(static_cast<unsigned int>(-1)) {}
  unsigned int index;
  Node *back;
  Node *next;
  unsigned int spid;
};


/**
 *  A structure for building a tree via stepwise addition.
 *  Provides methods to insert a taxon, to save the tree, and
 *  to access the tree topology.
 *  
 *  InducedStepwiseTree objects can be attached to the StepwiseTree
 *  using addListener. If so, those objects will be updated after
 *  each call to addLeaf
 *  
 */
class StepwiseTree
{
  public:
    StepwiseTree();

    virtual ~StepwiseTree();
 
    /**
     *  Insert a new leaf with the given spid into branch
     */
    Node *addLeaf(unsigned int spid, 
        Node *branch);

    /**
     *  Randomly return a branch
     */
    Node *getRandomBranch();

    /**
     *  Return a branch
     */ 
    Node *getAnyBranch() const;

    /**
     *  Return a leaf
     */
    Node *getAnyLeaf() const;

    /**
     *  Connects an InducedStepwiseTree to this
     */
    void addListener(InducedStepwiseTree *listener);

    /**
     *  Return the newick string corresponding to the tree 
     */
    std::string getNewickString(const std::vector<std::string> &spidToLabel) const;

    /**
     *  Return all the branches (if node is in the branches, node->back
     *  won't be in the branches)
     */
    std::vector<Node *> getBranches() const;

    /**
     *  Assuming that n2 corresponds to n1->next or n1->next->next,
     *  return the next node that is neither n1 nor n2
     */
    static Node *getOtherNext(Node *n1, Node *n2);
    
    /**
     *  Returns the vector of nodes of the current tree
     */
    const std::vector<Node *> &getAllNodes() const {return _nodes;}
  
    /**
     *  Return a vector that contains, for each listening induced tree,
     *  the induced node corresponding to the last inserted supernode.
     *  If the corresponding induced tree does not cover the last
     *  inserted supernode, the entry is set to nullptr
     */
    const std::vector<Node *> &getLastAddedInducedNodes() const {
      return _lastAddedInducedNode;
    }

    /**
     *  Fill nodes with all nodes under root (root is also inserted)
     *  in a post order traversal
     */
    static void fillPostOrder(Node *root, std::vector<Node *> &nodes);
  private:
    void addNode(Node *node);
    std::vector<Node *> _nodes;
    std::vector<InducedStepwiseTree*> _listeners;
    std::vector<Node *> _lastAddedInducedNode;
};

