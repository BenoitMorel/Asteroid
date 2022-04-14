#include "StepwiseTree.hpp"
#include <sstream>
#include <iostream>
#include <maths/Random.hpp>
#include <cassert>
#include <trees/InducedStepwiseTree.hpp>

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

StepwiseTree::StepwiseTree()
{
}

Node *StepwiseTree::addLeaf(unsigned int spid,
    Node *branch)
{
  auto leaf = new Node();
  leaf->spid = spid;
  addNode(leaf);
  if (_nodes.size() == 1) { // first  leaf
  } else if (_nodes.size() == 2) { // second leaf
    linkBack(leaf, branch);
  } else {
    auto n1 = new Node();
    auto n2 = new Node();
    auto n3 = new Node();
    addNode(n1);
    addNode(n2);
    addNode(n3);
    auto b1 = branch;
    auto b2 = branch->back;
    linkBack(leaf, n1);
    linkBack(n2, b1);
    linkBack(n3, b2);
    linkNext(n1, n2, n3);
  }
  for (auto listener: _listeners) {
    listener->addLeaf(leaf, spid);
  }
  return leaf;
}


StepwiseTree::~StepwiseTree()
{
  for (auto node: _nodes) {
    delete node;
  }
}


static void printAux(const std::vector<std::string> &spidToLabel,
    Node *node,
    std::stringstream &ss)
{
  if (node->next) {
    ss << "(";
    printAux(spidToLabel, node->next->back, ss);
    ss << ",";
    printAux(spidToLabel, node->next->next->back, ss);
    ss << ")";
    ss << ":" << node->index << "-" <<  node->next->index << "-" <<  node->next->next->index;
  } else {
    if (node->spid < spidToLabel.size()) {
      ss << spidToLabel[node->spid];
    } else {
      ss << "null";
    }
    
    ss << ":" << node->index;
  }
}

std::string StepwiseTree::getNewickString(const std::vector<std::string> &spidToLabel) const
{
  if (!_nodes.size()) {
    return std::string();
  } else if (_nodes.size() == 1) {
    std::string res;
    if (spidToLabel.size()) {
      res = spidToLabel.at(_nodes[0]->spid);
    }
    res += ";";
    return res;
  }
  std::stringstream ss;
  auto root = _nodes[1];
  ss << "(";
  printAux(spidToLabel, root, ss);
  ss << ",";
  printAux(spidToLabel, root->back, ss);
  ss << ");";
  return ss.str();
}
 

Node *StepwiseTree::getRandomBranch()
{
  if (!_nodes.size()) {
    return nullptr;
  }
  return _nodes[Random::getInt(_nodes.size())];
}
    
Node *StepwiseTree::getAnyBranch()
{
  if (!_nodes.size()) {
    return nullptr;
  } else {
    return _nodes[0];
  }
}
    
std::vector<Node *> StepwiseTree::getBranches() const
{
  std::vector<Node *> branches;
  if (_nodes.size() <= 1) {
    if (_nodes.size() == 1) {
      branches.push_back(_nodes[0]);
    }
    return branches;
  }
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
    
void StepwiseTree::addListener(InducedStepwiseTree *listener)
{
  _listeners.push_back(listener);
}

