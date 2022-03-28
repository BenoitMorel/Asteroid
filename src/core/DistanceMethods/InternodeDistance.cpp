#include "InternodeDistance.hpp"
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>

void computeFromGeneTreeAux(corax_unode_t *node,
    const IDParam &params,
    double d,
    std::vector<double> &distances)
{
  if (!node->next) {
    d += params.useBL ? node->length : 1.0;
    auto spid = static_cast<size_t>(
        reinterpret_cast<intptr_t>(node->data));
    distances[spid] = std::min(d, distances[spid]); 
    return;
  }
  bool isPolytomy = false;
  isPolytomy |= node->length  <= params.minBL;
  isPolytomy &= (nullptr != node->back->next);
  if (!isPolytomy) {
    d += params.useBL ? node->length : 1.0;
  }
  for (auto n = node->next; n != node; n = n->next) {
    computeFromGeneTreeAux(n->back, params, d, distances);
  }
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


void computeFromSpeciesTreeAux(corax_unode_t *node,
    double d,
    std::vector<double> &distances)
{
  d += 1.0;
  if (!node->next) {
    distances[node->node_index] = d; 
    return;
  }
  for (auto n = node->next; n != node; n = n->next) {
    computeFromSpeciesTreeAux(n->back, d, distances);
  }
}

void InternodeDistance::computeFromSpeciesTree(PLLUnrootedTree &speciesTree,
    DistanceMatrix &distanceMatrix)
{
  for (auto leaf: speciesTree.getLeaves()) {
    computeFromSpeciesTreeAux(leaf->back,
        0.0,
        distanceMatrix[leaf->node_index]);
  }
}



