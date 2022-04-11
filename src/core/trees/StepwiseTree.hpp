#pragma once

#include <string>
#include <vector>

struct Node {
  Node(): index(0), back(nullptr), next(nullptr), spid(0) {}
  unsigned int index;
  Node *back;
  Node *next;
  unsigned int spid;
};

class StepwiseTree
{
  public:
    StepwiseTree();

    ~StepwiseTree();
 
    void addLeaf(unsigned int spid, 
        Node *branch);

    Node *getRandomBranch();
    Node *getAnyBranch();

    std::string getNewickString(const std::vector<std::string> &spidToLabel);

    std::vector<Node *> getBranches();
  private:
    void addNode(Node *node);
    std::vector<Node *> _nodes;
};

