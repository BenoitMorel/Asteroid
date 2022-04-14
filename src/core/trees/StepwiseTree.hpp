#pragma once

#include <string>
#include <vector>

class InducedStepwiseTree;

struct Node {
  Node(): index(-1), back(nullptr), next(nullptr), spid(-1) {}
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
    Node *getAnyBranch();
    void addListener(InducedStepwiseTree *listener);

    std::string getNewickString(const std::vector<std::string> &spidToLabel) const;

    std::vector<Node *> getBranches() const;
    
    const std::vector<Node *> &getAllNodes() const {return _nodes;}
  private:
    void addNode(Node *node);
    std::vector<Node *> _nodes;
    std::vector<InducedStepwiseTree*> _listeners;
};

