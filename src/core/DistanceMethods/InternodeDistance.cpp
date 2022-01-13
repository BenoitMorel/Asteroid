#include "InternodeDistance.hpp"
#include <trees/PLLUnrootedTree.hpp>


void computeFromGeneTreeAux(corax_unode_t *node,
    const IDParam &params,
    double d,
    std::vector<double> &distances)
{
  if (!node->next) {
    d += 1.0;
    auto spid = static_cast<size_t>(
        reinterpret_cast<intptr_t>(node->data));
    distances[spid] = std::min(d, distances[spid]); 
    return;
  }
  if (node->length >= params.minBL // contract min branch length
      || !node->back->next) { // always add one at first level of recursion
    d += 1.0;
  }
  computeFromGeneTreeAux(node->next->back, params, d, distances);
  computeFromGeneTreeAux(node->next->next->back, params, d, distances);
}

void InternodeDistance::computeFromGeneTree(PLLUnrootedTree &geneTree,
    DistanceMatrix &distanceMatrix,
    const IDParam &params)
{
  for (auto leaf: geneTree.getLeaves()) {
    auto spid = static_cast<size_t>(
        reinterpret_cast<intptr_t>(leaf->data));
    computeFromGeneTreeAux(leaf->back,
        params,
        0.0,
        distanceMatrix[spid]);
  }
}



