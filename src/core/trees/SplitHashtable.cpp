#include "SplitHashtable.hpp"
#include <trees/PLLUnrootedTree.hpp>  
#include <corax/corax.h>
#include <stdio.h>
#include <stdlib.h> 
#include <IO/Logger.hpp>

static void insert(std::unordered_map<BitVector, unsigned int> &splits,
    const BitVector &split)
{
  auto it = splits.find(split);
  if (it == splits.end()) {
    splits.insert({split, 1});
  } else {
    it->second++;
  }
}

static BitVector addSplitFrom(corax_unode_t *node,
    const BitVector &emptySplit,
    std::unordered_map<BitVector, unsigned int> &splits)
{
  auto split = emptySplit;
  if (!node->next) {
    split.set(node->node_index);
  } else {
    auto n1 = node->next->back;
    auto n2 = node->next->next->back;
    split |= addSplitFrom(n1, emptySplit, splits);
    split |= addSplitFrom(n2, emptySplit, splits);
    insert(splits, split);
    insert(splits, ~split);
  }
  return split;
}


void SplitHashtable::addTree(PLLUnrootedTree &tree)
{
  BitVector emptySplit(tree.getLeavesNumber(), false);
  // start recursion from a leaf
  addSplitFrom(tree.getNode(0)->back,
      emptySplit,
      _splits); 
}

static char* intToString(int num) { 
  int size = log10(num) + 1; 
  if (num == 0) {
    size = 2;
  }
  auto str = static_cast<char*>(malloc(size*sizeof(char)));
  sprintf(str, "%d",num); 
  return str;
} 

static void writeSupport(corax_unode_t *node,
    int support)
{
  free(node->label);
  node->label = intToString(support);
}

static BitVector setSupportValue(corax_unode_t *node,
    const BitVector &emptySplit,
    const std::unordered_map<BitVector, unsigned int> &splits)
{
  auto split = emptySplit;
  if (!node->next) {
    split.set(node->node_index);
  } else {
    auto n1 = node->next->back;
    auto n2 = node->next->next->back;
    split |= setSupportValue(n1, emptySplit, splits);
    split |= setSupportValue(n2, emptySplit, splits);
    if (!node->back->next) { // first recursion call 
      return split;
    }
    int support = 0;
    auto it = splits.find(split);
    if (it != splits.end()) {
      support = it->second;
    }
    writeSupport(node, support);
    writeSupport(node->back, support);
  }
  return split;
}
  
void SplitHashtable::computeSupportValues(PLLUnrootedTree &tree)
{
  BitVector emptySplit(tree.getLeavesNumber(), false);
  // start the recursion from a leaf
  setSupportValue(tree.getNode(0)->back,
      emptySplit,
      _splits);
}

