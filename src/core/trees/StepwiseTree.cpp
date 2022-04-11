#include "StepwiseTree.hpp"
#include <sstream>
#include <iostream>
#include <maths/Random.hpp>
#include <cassert>

static void linkBack(Node *n1, Node *n2)
{
  n1->back = n2;
  n2->back = n1;
}

static void linkNext(Node *n1, Node *n2, Node *n3)
{
  n1->next = n2;
  n2->next = n3;
  n3->next = n1;
}

StepwiseTree::StepwiseTree(const std::string &label1,
    const std::string &label2,
    const std::string &label3,
    unsigned int spid1,
    unsigned int spid2,
    unsigned int spid3)
{
  for (unsigned int i = 0; i < 6; ++i) {
    addNode(new Node());
  }
  _nodes[0]->label = label1;
  _nodes[1]->label = label2;
  _nodes[2]->label = label3;
  _nodes[0]->spid = spid1;
  _nodes[1]->spid = spid2;
  _nodes[2]->spid = spid3;
  for (unsigned int i = 0; i < 3; ++i) {
    auto leaf = _nodes[i];
    auto node = _nodes[i + 3];
    linkBack(leaf, node);
    node->next = _nodes[3 + ((i+1) % 3)];
  }
}

StepwiseTree::~StepwiseTree()
{
  for (auto node: _nodes) {
    delete node;
  }
}


static void printAux(Node *node,
    std::stringstream &ss)
{
  if (node->next) {
    ss << "(";
    printAux(node->next->back, ss);
    ss << ",";
    printAux(node->next->next->back, ss);
    ss << ")";
  }
  ss << node->label << ":1.0";
}

std::string StepwiseTree::getNewickString()
{
  std::stringstream ss;
  auto root = _nodes[1];
  ss << "(";
  printAux(root, ss);
  ss << ",";
  printAux(root->back, ss);
  ss << ");";
  return ss.str();
}
 

void StepwiseTree::addLeaf(const std::string &label, 
    Node *branch,
    unsigned int spid)
{
  auto b1 = branch;
  auto b2 = branch->back;
  auto n1 = new Node();
  auto n2 = new Node();
  auto n3 = new Node();
  auto leaf = new Node();
  addNode(leaf);
  addNode(n1);
  addNode(n2);
  addNode(n3);
  leaf->label = label;
  leaf->spid = spid; 
  linkBack(leaf, n1);
  linkBack(n2, b1);
  linkBack(n3, b2);
  linkNext(n1, n2, n3);
}

Node *StepwiseTree::getRandomBranch()
{
  return _nodes[Random::getInt(_nodes.size())];
}
    
std::vector<Node *> StepwiseTree::getBranches()
{
  std::vector<Node *> branches;
  std::vector<bool> marked(_nodes.size(), false);
  for (auto branch: _nodes) {
    if (marked[branch->index]) {
      continue;
    }
    branches.push_back(branch);
    marked[branch->index] = true;
    marked[branch->back->index] = true;
  }
  assert(branches.size() == _nodes.size() / 2);
  return branches;
}

void StepwiseTree::addNode(Node *node)
{
  node->index = _nodes.size();
  _nodes.push_back(node);
}

