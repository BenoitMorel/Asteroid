#pragma once

#include <maths/bitvector.hpp>
#include <unordered_map>

using BitVector = genesis::utils::Bitvector; 
class PLLUnrootedTree;

class SplitHashtable {
public:
  SplitHashtable() {}
  virtual ~SplitHashtable() {}


  /**
   * Add splits from the input tree.
   * IMPORTANT: The splits are based on the 
   * tree leaf node indices
   */
  void addTree(PLLUnrootedTree &tree);

  /**
   * Compute the support values of the internal branches
   * of tree and stores them in the label field.
   * IMPORTANT: The splits are based on the 
   * tree leaf node indices
   */

  void computeSupportValues(PLLUnrootedTree &tree);
private:
  std::unordered_map<BitVector, unsigned int> _splits;
  

};
