#pragma once

#include <string>
#include <vector>

struct Node {
  Node(): back(nullptr), next(nullptr) {}
  unsigned int index;
  Node *back;
  Node *next;
  std::string label;
};

class StepwiseTree
{
  public:
    StepwiseTree(const std::string &label1,
        const std::string &label2,
        const std::string &label3);

    ~StepwiseTree();
 
    void addLeaf(const std::string &label, 
        Node *branch);

    Node *getRandomBranch();

    std::string getNewickString();

    std::vector<Node *> getBranches();
  private:
    void addNode(Node *node);
    std::vector<Node *> _nodes;
};

