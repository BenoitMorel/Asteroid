#pragma once

#include <util/types.hpp>

class PLLUnrootedTree;

struct IDParam {
  IDParam():minBL(-1.0) {}

  double minBL; // contract internal branches shorter than minBL
};

class InternodeDistance {
public:
  /**
   *  Compute the MiniNJ distance matrix between all species
   *  in the gene tree. The data field of the nodes should
   *  countain the spid of the corresponding species
   *
   *  distanceMatrix should already be allocated and filled
   *  with +inf
   */
  static void computeFromGeneTree(PLLUnrootedTree &geneTree,
    DistanceMatrix &distanceMatrix,
    const IDParam &param);

  static void computeFromSpeciesTree(PLLUnrootedTree &geneTree,
    DistanceMatrix &distanceMatrix);

};



