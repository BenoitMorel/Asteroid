#include "PLLUnrootedTree.hpp"

#include <IO/Logger.hpp>
#include <IO/LibpllException.hpp>
#include <stack>
#include <functional>
#include <sstream>
#include <deque>
#include <algorithm>
#include <map>

void defaultUnodePrinter(corax_unode_t *node, 
    std::stringstream &ss)
{
  if (node->label) {
    ss << node->label;
  }
  ss << ":" << node->length;
}


static void destroyNodeData(void *)
{
}

static void utreeDestroy(corax_utree_t *utree) {
  if(!utree)
    return;
  corax_utree_destroy(utree, destroyNodeData);
}



static corax_utree_t *readNewickFromStr(const std::string &str) 
{
  auto utree =  corax_utree_parse_newick_string_unroot(str.c_str());
  if (!utree) 
    throw LibpllException("Error while reading tree from std::string: ", str);
  return utree;
}


static corax_utree_t *readNewickFromFile(const std::string &str)
{
  try {
    auto utree =  corax_utree_parse_newick_unroot(str.c_str());
    if (!utree) 
      throw LibpllException("Error while reading tree from file: ", str);
    return utree;
  } catch (...) {
      throw LibpllException("Error while reading tree from file: ", str);
  }
}

static corax_utree_t *buildUtree(const std::string &str, bool isFile)
{
  if (isFile) {
    return readNewickFromFile(str);
  } else {
    return readNewickFromStr(str);
  }
}

static void setNodeLabel(corax_unode_t *node, const std::string &label)
{
  free(node->label);
  node->label = static_cast<char*>(malloc(sizeof(char) * (label.size() + 1)));
  std::strcpy(node->label, label.c_str());
}

void PLLUnrootedTree::setLabel(unsigned int nodeIndex, const std::string &label)
{
  setNodeLabel(getNode(nodeIndex), label);
}
  


PLLUnrootedTree::PLLUnrootedTree(const std::string &str, bool isFile):
  _tree(buildUtree(str, isFile), utreeDestroy)
{
}

PLLUnrootedTree::PLLUnrootedTree(const std::vector<const char*> &labels,
    unsigned int seed):
  _tree(corax_utree_random_create(static_cast<unsigned int>(labels.size()), &labels[0], seed), utreeDestroy)
{

}

static corax_utree_t *createTreeFromLabels(const std::unordered_set<std::string> &labels, 
    unsigned int seed)
{
  std::vector<const char *> temp;
  for (auto &label: labels) {
    temp.push_back(label.c_str());
  }
  return corax_utree_random_create(static_cast<unsigned int>(temp.size()),
        &temp[0], 
        seed);

}

PLLUnrootedTree::PLLUnrootedTree(const std::unordered_set<std::string> &labels,
    unsigned int seed):
  _tree(createTreeFromLabels(labels, seed), utreeDestroy)
{
}

void PLLUnrootedTree::save(const std::string &fileName)
{
  std::ofstream os(fileName, std::ofstream::out);
  char *newick = corax_utree_export_newick_rooted(getRawPtr()->nodes[0], 0);
  os << newick;
  os.close();
  free(newick);
}

void PLLUnrootedTree::setMissingBranchLengths(double minBL)
{
  for (auto node: getLeaves()) {
    if (0.0 == node->length) {
      node->length = minBL;
    } 
  }
  for (unsigned int i = _tree->tip_count; i < _tree->tip_count + _tree->inner_count; ++i) {
    if (0.0 == _tree->nodes[i]->length)
      _tree->nodes[i]->length = minBL;
    if (0.0 == _tree->nodes[i]->next->length)
      _tree->nodes[i]->next->length = minBL;
    if (0.0 == _tree->nodes[i]->next->next->length)
      _tree->nodes[i]->next->next->length = minBL;
  }  
}
  
void PLLUnrootedTree::setAllBranchLengths(double val)
{
  for (auto node: getLeaves()) {
    node->length = val;
  }
  for (unsigned int i = _tree->tip_count; i < _tree->tip_count + _tree->inner_count; ++i) {
    _tree->nodes[i]->length = val;
    _tree->nodes[i]->next->length = val;
    _tree->nodes[i]->next->next->length = val;
  }  
}

  
CArrayRange<corax_unode_t*> PLLUnrootedTree::getLeaves() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<corax_unode_t*> PLLUnrootedTree::getNodes() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes, getNodesNumber());
}

CArrayRange<corax_unode_t*> PLLUnrootedTree::getInnerNodes() const
{
  return CArrayRange<corax_unode_t*>(_tree->nodes + getLeavesNumber(), getInnerNodesNumber());
}


unsigned int PLLUnrootedTree::getNodesNumber() const
{
  return getLeavesNumber() + getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getDirectedNodesNumber() const
{
  return getLeavesNumber() + 3 * getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getLeavesNumber() const
{
  return _tree->tip_count;
}

unsigned int PLLUnrootedTree::getInnerNodesNumber() const
{
  return _tree->inner_count;
}
  
corax_unode_t *PLLUnrootedTree::getNode(unsigned int node_index) const
{
  return _tree->nodes[node_index];
}

corax_unode_t *PLLUnrootedTree::getAnyInnerNode() const
{
  return getNode(getLeavesNumber());
}
  
std::unordered_set<std::string> PLLUnrootedTree::getLeavesLabels()
{
  std::unordered_set<std::string> res;
  for (auto leaf: getLeaves()) {
    if (leaf->label) {
      res.insert(std::string(leaf->label));
    }
  }
  return res;
}


static bool isBranchIn(corax_unode_t *b, 
    const std::unordered_set<corax_unode_t *> &branches)
{
  return branches.find(b) != branches.end() 
    || branches.find(b->back) != branches.end();
}

std::unordered_set<corax_unode_t *> PLLUnrootedTree::getBranches() const
{
  std::unordered_set<corax_unode_t *> branches;
  for (auto node: getNodes()) {
    if (!isBranchIn(node, branches)) {
      branches.insert(node);
    }
    if (node->next) {
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
    }
  }
  return branches;
}

struct ToSort {
  corax_unode_t *node;
  std::string label;
  unsigned int distanceToLeaf;
 
  ToSort() {}
  ToSort(corax_unode_t *inode, 
      const std::string &ilabel,
      unsigned int idistanceToLeaf): node(inode),
  label(ilabel),
  distanceToLeaf(idistanceToLeaf) {}
  ToSort(const ToSort &s): node(s.node),
    label(s.label),
    distanceToLeaf(s.distanceToLeaf) {}
};

static void fillBranchesRec(corax_unode_t * node,
    std::vector<ToSort> &toSortVector,
    ToSort &toSort)
{
  if (!node->next) {
    toSort = ToSort(node, std::string(node->label), 0);
    toSortVector.push_back(toSort);
    return;
  }
  ToSort toSort1;
  fillBranchesRec(node->next->back, toSortVector, toSort1);
  ToSort toSort2;
  fillBranchesRec(node->next->next->back, toSortVector, toSort2);
  if (toSort1.label > toSort2.label) {
    toSort = ToSort(node, toSort1.label, toSort1.distanceToLeaf + 1);
  } else {
    toSort = ToSort(node, toSort2.label, toSort2.distanceToLeaf + 1);
  }
  toSortVector.push_back(toSort);
}

struct less_than_key
{
  inline bool operator() (const ToSort& t1, const ToSort& t2)
  {
    if (t1.label == t2.label) {
      return t1.distanceToLeaf < t2.distanceToLeaf;
    } 
    return t1.label < t2.label;
  }
};


static void fillPostOrder(corax_unode_t *node,
    std::vector<corax_unode_t*> &nodes,
    std::vector<bool> &markedNodes)
{
  // we already traversed this node
  if (markedNodes[node->node_index]) {
    return;
  }
  // mark the node as traversed
  markedNodes[node->node_index] = true;
  // first process children
  if (node->next) {
    fillPostOrder(node->next->back, nodes, markedNodes);
    fillPostOrder(node->next->next->back, nodes, markedNodes);
  }
  nodes.push_back(node);
}

static void fillPostOrderNoMarker(corax_unode_t *node,
    std::vector<corax_unode_t*> &nodes)
{

  if (node->next) {
    fillPostOrderNoMarker(PLLUnrootedTree::getLeft(node), nodes);
    fillPostOrderNoMarker(PLLUnrootedTree::getRight(node), nodes);
  }
  nodes.push_back(node);
}
  
std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodesFrom(corax_unode_t *node) const
{
  std::vector<corax_unode_t*> nodes;
  fillPostOrderNoMarker(node, nodes);
  return nodes;
}


std::vector<corax_unode_t*> PLLUnrootedTree::getPostOrderNodes(bool innerOnly) const
{
  std::vector<corax_unode_t*> nodes;
  std::vector<bool> markedNodes(getDirectedNodesNumber(), false);
  if (innerOnly) {
    for (auto node: getLeaves()) {
      markedNodes[node->node_index] = true;
    }
  }
  // do the post order traversal from all possible virtual roots 
  for (auto node: getLeaves()) {
    fillPostOrder(node->back, nodes, markedNodes);
  }
  if (innerOnly) {
    assert(nodes.size() == getDirectedNodesNumber() - getLeavesNumber());
  } else {
    assert(nodes.size() == getDirectedNodesNumber());
  }
  return nodes;
}
  
std::vector<corax_unode_t*> PLLUnrootedTree::getReverseDepthNodes() const
{
  std::deque<corax_unode_t *> q;
  std::vector<bool> marked(getDirectedNodesNumber(), false);
  std::vector<corax_unode_t *> nodes;
  for (auto leaf: getLeaves()) {
    marked[leaf->node_index] = true;
    q.push_back(leaf);
  }
  while (!q.empty()) {
    auto node = q.front();
    q.pop_front();
    nodes.push_back(node);
    auto back = node->back;
    if (!back->next) {
      continue;
    }
    auto left = back->next;
    auto right = back->next->next;
    if (!marked[left->node_index]) {
      q.push_back(left);
      marked[left->node_index] = true;
    }
    if (!marked[right->node_index]) {
      q.push_back(right);
      marked[right->node_index] = true;
    }
  }
  assert(nodes.size() == getDirectedNodesNumber());
  return nodes; 
}

static void computePairwiseDistancesRec(corax_unode_t *currentNode, 
    double currentDistance,
    VectorDouble &distancesVector)
{
  currentDistance += currentNode->length;
  if (!currentNode->next) {
    // leaf
    distancesVector[currentNode->node_index] = currentDistance;
    return;
  }
  computePairwiseDistancesRec(currentNode->next->back, 
      currentDistance, 
      distancesVector);
  computePairwiseDistancesRec(currentNode->next->next->back, 
      currentDistance, 
      distancesVector);
} 

void PLLUnrootedTree::computePairwiseDistances(MatrixDouble &distances,
    bool leavesOnly)
{
  auto M = leavesOnly ? getLeavesNumber() : getDirectedNodesNumber(); 
  auto N = getLeavesNumber();
  auto MIter = leavesOnly ? getLeavesNumber() : getNodesNumber();
  VectorDouble zeros(N, 0.0);
  distances = MatrixDouble(M, zeros);
  for (unsigned int i = 0; i < MIter; ++i) {
    auto node = getNode(i);
    auto &distancesVector = distances[node->node_index];
    computePairwiseDistancesRec(node->back, 0.0, distancesVector);
    if (node->next) {
      // compute distances to leaves in all three directions
      computePairwiseDistancesRec(node->next->back, 
          0.0, 
          distancesVector);
      computePairwiseDistancesRec(node->next->next->back, 
          0.0, 
          distancesVector);
      // also update the two other directed nodes
      distances[node->next->node_index] = distancesVector;
      distances[node->next->next->node_index] = distancesVector;
    }
  }
}


static void getCladeRec(corax_unode_t *node, 
    std::unordered_set<unsigned int> &clade)
{
  if (node->next) {
    getCladeRec(node->next->back, clade);
    getCladeRec(node->next->next->back, clade);
  } else {
    clade.insert(node->node_index);
  }
}

std::unordered_set<unsigned int> 
  PLLUnrootedTree::getClade(corax_unode_t *node)
{
  std::unordered_set<unsigned int> clade;
  getCladeRec(node, clade);
  return clade;
}




std::vector<double> PLLUnrootedTree::getMADRelativeDeviations()
{
  MatrixDouble distances;
  computePairwiseDistances(distances, false);
  auto nodes = getPostOrderNodes();  
  std::vector<double> deviations(nodes.size());
  for (auto node: getPostOrderNodes()) {
    // we reuse notations from the MAD paper (I, J, i, j, rho, b, c)
    // I and J are the set of leaves of each side of the branch node
    // rho is the relative position of the rooting that minimizes 
    // the squared relative deviation
    auto I = getClade(node); 
    auto J = getClade(node->back);
    auto rho = 0.0;
    auto rhoDen = 0.0;
    auto i = node->node_index;
    auto Dij = node->length;
    for (auto b: I) {
      for (auto c: J) {
        auto invPowBC = pow(distances[b][c], -2.0);
        rho += (distances[b][c] - 2.0 * distances[i][b]) * invPowBC;
        rhoDen += invPowBC;
      }
    }
    rho = rho / (2.0 * Dij * rhoDen);
    rho = std::max(std::min(rho, 1.0), 0.0);
    auto deviation = 0.0;
    for (auto b: I) {
      for (auto c: J) {
        auto v = (2.0 * (distances[i][b] + rho * Dij) / distances[b][c]);
        v -= 1.0;
        deviation += pow(v, 2.0);
      }
    }
    deviations[node->node_index] = deviation;
  }
  return deviations;
}

static void printAux(corax_unode_t *node,
    std::stringstream &ss,
    UnodePrinter f)
{
  if (node->next) {
    ss << "(";
    printAux(node->next->back, ss, f);
    ss << ",";
    printAux(node->next->next->back, ss, f);
    ss << ")";
  }
  f(node, ss);
}
  
std::string PLLUnrootedTree::getSubtreeString(corax_unode_t *subtree, UnodePrinter f)
{
  if (!subtree) {
    return "null";
  }
  std::stringstream ss;
  printAux(subtree, ss, f);
  return ss.str();
}

std::string PLLUnrootedTree::getNewickString(UnodePrinter f,
      corax_unode_t *root, 
      bool rooted) const
{
  std::stringstream ss;
  if (!root) {
    root = getAnyInnerNode();
  }
  if (rooted) {
    ss << "(";
    printAux(root, ss, f);
    ss << ",";
    printAux(root->back, ss, f);
    ss << ");";
  } else {
    ss << "(";
    printAux(root->back, ss, f);
    ss << ",";
    printAux(root->next->back, ss, f);
    ss << ",";
    printAux(root->next->next->back, ss, f);
    ss << ");";
  }
  return ss.str();
}

std::unordered_set<std::string> PLLUnrootedTree::getLabels() const
{
  std::unordered_set<std::string> res;
  for (auto node:  getLeaves()) {
    if (node->label) {
      res.insert(node->label);
    }
  }
  return res;
}

static corax_unode_t *searchForSet(corax_unode_t *node, 
    std::unordered_set<std::string> &currentNodeSet,
    const std::unordered_set<std::string> &set)
{
  if (!node->next) {
    if (node->label) {
      currentNodeSet.insert(std::string(node->label));
    }
  } else {
    auto left = node->next->back;
    auto right = node->next->next->back;
    std::unordered_set<std::string> rightNodeSet;
    auto res1 = searchForSet(left, currentNodeSet, set);
    if (res1) {
      return res1;
    }
    auto res2 = searchForSet(right, rightNodeSet, set);
    if (res2) {
      return res2;
    }
    for (auto &elem: rightNodeSet) {
      currentNodeSet.insert(elem);
    }
  }
  return (currentNodeSet == set) ? node : nullptr;
}



static size_t leafHash(corax_unode_t *leaf) {
  assert(leaf);
  std::hash<std::string> hash_fn;
  return hash_fn(std::string(leaf->label));
}

static size_t getTreeHashRec(corax_unode_t *node, size_t i) {
  assert(node);
  if (i == 0) 
    i = 1;
  if (!node->next) {
    return leafHash(node);
  }
  auto hash1 = getTreeHashRec(node->next->back, i + 1);
  auto hash2 = getTreeHashRec(node->next->next->back, i + 1);
  std::hash<size_t> hash_fn;
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  return hash_fn(m * i + M);

}

static corax_unode_t *findMinimumHashLeafRec(corax_unode_t * root, size_t &hashValue)
{
  assert(root);
  if (!root->next) {
    hashValue = leafHash(root);
    return root;
  }
  auto n1 = root->next->back;
  auto n2 = root->next->next->back;
  size_t hash1, hash2;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    hashValue = hash1;
    return min1;
  } else {
    hashValue = hash2;
    return min2;
  }
}

static corax_unode_t *findMinimumHashLeaf(corax_unode_t * root) 
{
  assert(root);
  auto n1 = root;
  auto n2 = root->back;
  size_t hash1 = 0;
  size_t hash2 = 0;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    return min1;
  } else {
    return min2;
  }
}

size_t PLLUnrootedTree::getUnrootedTreeHash() const
{
  auto minHashLeaf = findMinimumHashLeaf(getAnyInnerNode());
  auto res = getTreeHashRec(minHashLeaf, 0) + getTreeHashRec(minHashLeaf->back, 0);
  return res;
}



static const char *orderChildren(corax_unode_t *node,
    std::vector<bool> &leftFirst)
{
  if (!node->next) {
    return node->label;
  }
  const char *label1 = orderChildren(node->next->back, leftFirst);
  const char *label2 = orderChildren(node->next->next->back, leftFirst);
  if (strcmp(label1, label2) < 0) {
    leftFirst[node->node_index] = true;
    return label1;
  } else {
    leftFirst[node->node_index] = false;
    return label2;
  }
}

static bool areIsomorphicAux(corax_unode_t *node1,
    corax_unode_t *node2,
    const std::vector<bool> &leftFirst1,
    const std::vector<bool> &leftFirst2)
{
  if (!node1->next || !node2->next) {
    // at least one is a leaf
    if (node1->next || node2->next) {
      // only one is a leaf 
      return false;
    }
    return strcmp(node1->label, node2->label) == 0;
  }
  // both are internal nodes
  auto l1 = node1->next->back;
  auto r1 = node1->next->next->back;
  auto l2 = node2->next->back;
  auto r2 = node2->next->next->back;
  if (!leftFirst1[node1->node_index]) {
    std::swap(l1, r1);
  }
  if (!leftFirst2[node2->node_index]) {
    std::swap(l2, r2);
  }
  return areIsomorphicAux(l1, l2, leftFirst1, leftFirst2) 
    && areIsomorphicAux(r1, r2, leftFirst1, leftFirst2);
}

bool PLLUnrootedTree::areIsomorphic(const PLLUnrootedTree &t1,
    const PLLUnrootedTree &t2)
{
  if (t1.getNodesNumber() != t2.getNodesNumber()) {
    return false;
  }
  auto startingLeaf1 = findMinimumHashLeaf(t1.getAnyInnerNode())->back;
  auto startingLeaf2 = findMinimumHashLeaf(t2.getAnyInnerNode())->back;
  std::vector<bool> leftFirst1(t1.getDirectedNodesNumber(), true);
  std::vector<bool> leftFirst2(t2.getDirectedNodesNumber(), true);
  orderChildren(startingLeaf1, leftFirst1);
  orderChildren(startingLeaf2, leftFirst2);
  return areIsomorphicAux(startingLeaf1, 
      startingLeaf2, 
      leftFirst1, 
      leftFirst2);
}

bool PLLUnrootedTree::isBinary() const
{
  for (auto node: getInnerNodes()) {
    assert(node->next);
    if (node->next->next->next != node) {
      return false;
    }
  }
  return true;
}
  
corax_unode_t *PLLUnrootedTree::findLeaf(const std::string &label)
{
  for (auto leaf: getLeaves()) {
    if (label == leaf->label) {
      return leaf;
    }
  }
  return nullptr;
}

void PLLUnrootedTree::ensureUniqueLabels() 
{
  auto labels = getLabels();
  unsigned int i = 0;
  std::string prefix("s");
  for (auto node: getPostOrderNodes()) {
    if (node->next) {
      std::string newLabel;
      if (node->label) {
        newLabel = std::string(node->label);
      }
      while (labels.find(newLabel) != labels.end() || newLabel.size() == 0) {
        newLabel = prefix + std::to_string(i++);
      }
      free(node->label);
      node->label = static_cast<char*>(malloc(sizeof(char) * (newLabel.size() + 1)));
      std::strcpy(node->label, newLabel.c_str());
      labels.insert(newLabel);
    }
  }
}
  
void PLLUnrootedTree::reindexLeaves(const StringToUint &labelToIndex)
{
  std::vector<corax_unode_t *> copy;
  for (auto leaf: getLeaves()) {
    copy.push_back(leaf);
  }
  for (auto leaf: copy) {
    auto spid = labelToIndex.at(leaf->label);
    leaf->node_index = spid;
    _tree->nodes[spid] = leaf;
  }
}
 
void PLLUnrootedTree::sortLeaves()
{
  std::map<std::string, corax_unode_t *> labelToNode;
  for (auto leaf: getLeaves()) {
    labelToNode.insert({std::string(leaf->label), leaf});
  }
  unsigned int i = 0;
  for (auto it: labelToNode) {
    auto node = it.second;
    _tree->nodes[i] = node;
    node->node_index = i;
    i++;
  }
}
 
static std::string getInducedNewickRec(corax_unode_t *node,
    const std::vector<bool> &coverage)
{
  if (!node->next) {
    if (coverage[node->node_index]) {
      return std::string(node->label);
    } else {
      return std::string();
    }
  }
  auto strLeft = getInducedNewickRec(node->next->back, coverage);
  auto strRight = getInducedNewickRec(node->next->next->back, coverage);
  if (!strLeft.size() && !strRight.size()) {
    return std::string();
  } 
  if (!strLeft.size()) {
    return strRight;
  }
  if (!strRight.size()) {
    return strLeft;
  }
  std::string res;
  res += "(";
  res += strLeft;
  res += ",";
  res += strRight;
  res += ")";
  return res;
}

std::string PLLUnrootedTree::getInducedNewick(const std::vector<bool> &coverage) const
{
  std::stringstream ss;
  unsigned int first = 0;
  for (; first < coverage.size(); ++first) {
    if (coverage[first]) {
      break;
    }
  }
  assert(first != coverage.size());
  auto firstLeaf = getNode(first);
  ss << "(";
  ss << firstLeaf->label;
  ss << ",";
  ss << getInducedNewickRec(firstLeaf->back, coverage);
  ss << ");";
  return ss.str();
}

std::shared_ptr<PLLUnrootedTree> PLLUnrootedTree::getInducedTree(const std::vector<bool> &coverage) const
{
  auto str = getInducedNewick(coverage);
  return std::make_shared<PLLUnrootedTree>(str, false);
}

static corax_unode_t *getOtherNext(corax_unode_t *n1, 
    corax_unode_t *n2)
{
  assert(n1 != n2);
  if (n1->next == n2) {
    assert(n2->next->next == n1);
    return n2->next;
  } else {
    assert(n1->next->next == n2);
    return n1->next;
  }
}
  
std::shared_ptr<PLLUnrootedTree> PLLUnrootedTree::getInducedTree(const BitVector &okNodeIndice) const
{
  std::vector<bool> v(okNodeIndice.size(), false);
  for (unsigned int i = 0; i < okNodeIndice.size(); ++i) {
    v[i] = okNodeIndice[i];
  }
  return getInducedTree(v);
}

static void fillWithChildren(corax_unode_t *super, 
    corax_unode_t *induced,
    std::vector<corax_unode_t *> &superToInducedRegraft)
{
  superToInducedRegraft[super->node_index] = induced;
  if (super->next) {
    fillWithChildren(PLLUnrootedTree::getLeft(super), induced, superToInducedRegraft);
    fillWithChildren(PLLUnrootedTree::getRight(super), induced, superToInducedRegraft);
  }
}

// has to be called in post order (wrt the super tree) traversal
static void mapNodesWithInducedTreeAux(corax_unode_t *superNode,
    const LabelToLeaf &inducedLabelToLeaf,  
    NodeVector &superToInduced,
    NodeVector &superToInducedRegraft)
{
  if (!superNode->next) {
    auto it = inducedLabelToLeaf.find(std::string(superNode->label));
    if (it != inducedLabelToLeaf.end()) {
      auto inducedLeaf = it->second;
      superToInduced[superNode->node_index] = inducedLeaf;
      superToInducedRegraft[superNode->node_index] = inducedLeaf;
    }
    return;
  }
  auto superLeft = PLLUnrootedTree::getLeft(superNode);
  auto superRight = PLLUnrootedTree::getRight(superNode);
  auto inducedLeft = superToInduced[superLeft->node_index]; 
  auto inducedRight = superToInduced[superRight->node_index]; 
  if (!inducedRight && !inducedLeft) {
    return;
  } else if (!inducedRight) {
    superToInduced[superNode->node_index] = inducedLeft;
    superToInducedRegraft[superNode->node_index] = inducedLeft;
    fillWithChildren(superRight, 
        inducedLeft,
        superToInducedRegraft);
  } else if (!inducedLeft) {
    superToInduced[superNode->node_index] = inducedRight;
    superToInducedRegraft[superNode->node_index] = inducedRight;
    fillWithChildren(superLeft, 
        inducedRight,
        superToInducedRegraft);
  } else {
    if (inducedLeft->back == inducedRight) {
      superToInducedRegraft[superNode->node_index] = inducedLeft;
      return;
    }
    auto inducedParent = getOtherNext(inducedLeft->back, inducedRight->back);
    superToInduced[superNode->node_index] = inducedParent;
    superToInducedRegraft[superNode->node_index] = inducedParent;
  }
}
    

void PLLUnrootedTree::mapNodesWithInducedTree(const PLLUnrootedTree &inducedTree,
      const NodeVector &superPostOrderNodes,
      NodeVector &superToInduced,
      NodeVector &superToInducedRegraft
      ) const
{  
  superToInduced.resize(superPostOrderNodes.size());
  superToInducedRegraft.resize(superPostOrderNodes.size());
  std::fill(superToInduced.begin(), superToInduced.end(), nullptr);
  std::fill(superToInducedRegraft.begin(), 
      superToInducedRegraft.end(),
      nullptr);
  LabelToLeaf inducedLabelToLeaf;
  for (auto leaf: inducedTree.getLeaves()) {
    inducedLabelToLeaf.insert({std::string(leaf->label), leaf});
  }
  for (auto superNode: superPostOrderNodes) {
    mapNodesWithInducedTreeAux(superNode,
        inducedLabelToLeaf,
        superToInduced,
        superToInducedRegraft);
  }
  for (auto superNode: superPostOrderNodes) {
    auto inducedRegraft = superToInducedRegraft[superNode->node_index];
    if (inducedRegraft) {
      superToInducedRegraft[superNode->back->node_index] = inducedRegraft;
    }
  }
}



