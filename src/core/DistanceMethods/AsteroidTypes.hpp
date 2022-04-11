#pragma once
#include <trees/PLLUnrootedTree.hpp>
#include <maths/bitvector.hpp>
#include <util/types.hpp>

using TreePtr = std::shared_ptr<PLLUnrootedTree>;
using Trees = std::vector<TreePtr>;
using GeneTrees = Trees;
using SpeciesTrees = Trees;
struct ScoredTree {
  TreePtr tree;
  double score;
  bool operator < (const ScoredTree& other) const {
    return other.score < score;
  }
};

using ScoredTrees = std::vector<ScoredTree>;

struct GeneCell {
  TreePtr geneTree;
  DistanceMatrix distanceMatrix;
  BitVector coverage;
  std::vector<unsigned int> gidToSpid;
  GeneCell(TreePtr ptr = nullptr): geneTree(ptr) {}
};


