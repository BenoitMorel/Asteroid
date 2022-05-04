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

class StepwiseTree
{
  public:
    StepwiseTree();

    virtual ~StepwiseTree();
 
    Node *addLeaf(unsigned int spid, 
        Node *branch);

    Node *getRandomBranch();
    Node *getAnyBranch() const;
    Node *getAnyLeaf() const;
    void addListener(InducedStepwiseTree *listener);

    std::string getNewickString(const std::vector<std::string> &spidToLabel) const;

    std::vector<Node *> getBranches() const;
    static Node *getOtherNext(Node *n1, Node *n2);
    const std::vector<Node *> &getAllNodes() const {return _nodes;}
  private:
    void addNode(Node *node);
    std::vector<Node *> _nodes;
    std::vector<InducedStepwiseTree*> _listeners;
};

