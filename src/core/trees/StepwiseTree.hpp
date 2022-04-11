#pragma once

#include <string>
#include <vector>

struct Node {
  Node(): index(0), back(nullptr), next(nullptr), spid(0) {}
  unsigned int index;
  Node *back;
  Node *next;
  std::string label;
  unsigned int spid;
};

class StepwiseTree
{
  public:
    StepwiseTree(const std::string &label1,
        const std::string &label2,
        const std::string &label3,
        unsigned int spid1 = -1,
        unsigned int spid2 = -1,
        unsigned int spid3 = -1);

    ~StepwiseTree();
 
    void addLeaf(const std::string &label, 
        Node *branch,
        unsigned int spid = -1);

    void removeLeaf(const std::string &label);

    Node *getRandomBranch();

    std::string getNewickString();

    std::vector<Node *> getBranches();
  private:
    void addNode(Node *node);
    std::vector<Node *> _nodes;
};

