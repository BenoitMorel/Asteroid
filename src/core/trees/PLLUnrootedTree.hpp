#pragma once

#include <corax/corax.h>
#include <unordered_set>
#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <util/types.hpp>
#include <util/CArrayRange.hpp>
#include <maths/Random.hpp>
#include <functional>
#include <cstring>
#include <maths/bitvector.hpp>


using UnodePrinter = 
  std::function<void(corax_unode_t *, std::stringstream &)>;
using NodeSet = std::unordered_set<corax_unode_t *>;
using NodeVector = std::vector<corax_unode_t *>;
using LabelToLeaf = std::unordered_map<std::string, corax_unode_t *>;

void defaultUnodePrinter(corax_unode_t *node, 
    std::stringstream &ss);

/**
 *  C++ wrapper around the libpll corax_utree_t structure
 *  to represent an unrooted tree
 */
class PLLUnrootedTree {
public:
  /**
   * Construct from a string that is either a path 
   * to a newick file or a newick string
   */
  explicit PLLUnrootedTree(const std::string &str, bool isFile = true);



  std::shared_ptr<PLLUnrootedTree> getInducedTree(const std::vector<bool> &okNodeIndices) const;
  std::shared_ptr<PLLUnrootedTree> getInducedTree(const BitVector &okNodeIndice) const;

  /**
   *  Construct a random tree from a set of taxa labels
   */
  explicit PLLUnrootedTree(const std::vector<const char*> &labels, unsigned int seed = static_cast<unsigned int>(Random::getInt()));
  explicit PLLUnrootedTree(const std::unordered_set<std::string> &labels,
      unsigned int seed = static_cast<unsigned int>(Random::getInt()));

  /**
   * Forbid copy
   */
  PLLUnrootedTree(const PLLUnrootedTree &) = delete;
  PLLUnrootedTree & operator = (const PLLUnrootedTree &) = delete;
  PLLUnrootedTree(PLLUnrootedTree &&) = delete;
  PLLUnrootedTree & operator = (PLLUnrootedTree &&) = delete;

  bool operator ==(const PLLUnrootedTree &other) const
  {
    return areIsomorphic(*this, other);
  }
  
  void setLabel(unsigned int nodeIndex, const std::string &label);
  
  /*
   * Tree dimension
   */
  // number of nodes in the node buffer (#leaves + #innner)
  unsigned int getNodesNumber() const;
  // #leaves + 3 * #inner
  unsigned int getDirectedNodesNumber() const;
  // #leaves
  unsigned int getLeavesNumber() const;
  // #inner
  unsigned int getInnerNodesNumber() const; 
  
  std::unordered_set<std::string> getLabels() const;

  /*
   * Node access
   */
  corax_unode_t *getAnyInnerNode() const;
  corax_unode_t *getNode(unsigned int node_index) const;

 
  
  /**
   * labels
   */
  std::unordered_set<std::string> getLeavesLabels();

  /*
   * Save the tree in newick format in filename
   */
  void save(const std::string &fileName); 
  std::string getNewickString(UnodePrinter f = defaultUnodePrinter,
      corax_unode_t *root = nullptr, 
      bool rooted = false) const;
  static std::string getSubtreeString(corax_unode_t *subtree, UnodePrinter f = defaultUnodePrinter);
  
  std::string getInducedNewick(const std::vector<bool> &coverage) const;

  /**
   * Replace null branch lengths with minBL
   */
  void setMissingBranchLengths(double minBL = 0.1); 
  
  /**
   *  Set all branch lengths to val
   */
  void setAllBranchLengths(double val = 1.0); 

  
  size_t getUnrootedTreeHash() const;

  /**
   *  Direct access to the libpll structure
   */
  corax_utree_t *getRawPtr() {return _tree.get();}
  
  CArrayRange<corax_unode_t*> getLeaves() const;

  /*
   *  C++11 range for accessing internal nodes. 
   *  Only includes one of the three pll elements per node
   *  in the tree
   */
  CArrayRange<corax_unode_t*> getInnerNodes() const;

  /*
   *  C++11 range for accessing nodes. 
   *  For internal nodes, Only includes one of the three pll elements per node
   *  in the tree
   */
  CArrayRange<corax_unode_t*> getNodes() const;

  /**
   *  Create a vector of all nodes (including all three pll
   *  internal elements per internal node), such that a node
   *  always comes after its (virtual) children.
   */
  std::vector<corax_unode_t*> getPostOrderNodes(bool innerOnly = false) const;
  std::vector<corax_unode_t*> getPostOrderNodesFrom(corax_unode_t *node) const;
  std::vector<corax_unode_t*> getReverseDepthNodes() const;
  
 
  /**
   *  Return a set of all branches (for each node, one and only 
   *  one of node and node->back will be inserted in the set)
   */
  std::unordered_set<corax_unode_t *> getBranches() const;

  /**
   *  Compute a matrix of pairwise distances.
   *  The distance between two nodes is the sum of the
   *  branch lengths along the path between these two nodes.
   *  
   *  distances[i][j] is the distance between node1 and node2,
   *  such that node1->node_index == i && node2->node_index == j
   *
   *  If leavesOnly is true, the output matrix is symetric, 
   *  and the diagonal elements are set to 0.0
   *
   *  If leavesOnly is false, distances is a rectangular matrix
   *  with dimensions distances[max_gene_id][leaves_number]
   *
   *  
   */
  void computePairwiseDistances(MatrixDouble &distances,
      bool leavesOnly = true);


  // TODO: DOCUMENT
  std::vector<double> getMADRelativeDeviations();

  /**
   *  Return the set of leaves under the input directed node
   */
  static std::unordered_set<unsigned int> getClade(corax_unode_t *node);

  static bool areIsomorphic(const PLLUnrootedTree &t1,
    const PLLUnrootedTree &t2);

  bool isBinary() const;
  
  void ensureUniqueLabels();

  /**
   *  Reorganize the node_index values of the leaves and
   *  reorder them in the utree::nodes buffer
   */
  void reindexLeaves(const StringToUint &labelToIndex);
  /**
   *  Sort the leaves by lexicographic order on their labels
   */
  void sortLeaves();

  /**
   *  Return the leaf node that has the label label,
   *  or null pointer if such leaf does not exist
   */
  corax_unode_t *findLeaf(const std::string &labe);

  static corax_unode_t *getLeft(corax_unode_t *node) {return node->next->back;}
  static corax_unode_t *getRight(corax_unode_t *node) {return node->next->next->back;}

  void mapNodesWithInducedTree(const PLLUnrootedTree &inducedTree,
      const NodeVector &superPostOrderNodes,
      NodeVector &superToInduced,
      NodeVector &superToInducedRegraft
      ) const;

private:
  std::unique_ptr<corax_utree_t, void(*)(corax_utree_t*)> _tree;
};





