#include <iostream>
#include <corax/corax.h>

#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <maths/Random.hpp>
#include "Arguments.hpp"
#include <trees/PLLUnrootedTree.hpp>
#include <trees/SplitHashtable.hpp>
#include <DistanceMethods/InternodeDistance.hpp>
#include <DistanceMethods/AsteroidOptimizer.hpp>
#include <memory>
#include <set>
#include <algorithm>
#include <IO/ParallelOfstream.hpp>
#include <maths/bitvector.hpp>
#include <trees/StepwiseTree.hpp>

using TreePtr = std::shared_ptr<PLLUnrootedTree>;
using Trees = std::vector<TreePtr>;
using GeneTrees = Trees;
using SpeciesTrees = Trees;
using Coverages = std::vector<BitVector>;

/**
 *  Read all trees from inputGeneTreeFile
 */
void parseTrees(const std::string &inputGeneTreeFile,
    GeneTrees &geneTrees)
{
  // Read all trees
  std::ifstream reader(inputGeneTreeFile);
  std::string geneTreeStr;
  while (std::getline(reader, geneTreeStr)) {
    if (geneTreeStr.size() < 2) {
      continue;
    }
    geneTrees.push_back(std::make_shared<PLLUnrootedTree>(
          geneTreeStr, false));
  }
}

 
/**
 *  Initialize the program state
 */
void init(Arguments &arg)
{
  ParallelContext::init(0);
  Logger::init();
  Logger::timed << "Starting Asteroid..." << std::endl;
  arg.printCommand();
  Random::setSeed(static_cast<unsigned int>(arg.seed));
}

/**
 *  Closes the program state
 */
void close()
{
  Logger::close();
  ParallelContext::finalize();
}

void readGeneTrees(const std::string &inputGeneTreeFile,
    GeneTrees &geneTrees)
{
  Logger::timed << "Parsing gene trees..." << std::endl;
  parseTrees(inputGeneTreeFile,
      geneTrees);
}

void getPerCoreGeneTrees(const Trees &geneTrees,
    Trees &perCoreGeneTrees)
{
  auto begin = ParallelContext::getBegin(geneTrees.size());
  auto end = ParallelContext::getEnd(geneTrees.size());
  for (unsigned int k = begin; k < end; ++k) {
    perCoreGeneTrees.push_back(geneTrees[k]);
  }
}

/**
 *  Fills the gene-to-species mapping object:
 *  - from the input mapping file
 *  - or from the gene tree leaf labels
 *  Identifies all the species labels and assign
 *  them a SPID
 *
 */
void extractMappings(const Arguments &arg,
    const GeneTrees &geneTrees,
    GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid)
{
  Logger::timed << "Mapping gene to species..." << std::endl;
  if (!arg.geneSpeciesMapping.size()) {
    for (auto geneTree: geneTrees) {
      mappings.fillFromGeneTree(*geneTree);
    }
  } else {
    mappings.fillFromMappingFile(arg.geneSpeciesMapping); 
  }
  std::set<std::string> sortedSpecies(mappings.getCoveredSpecies().begin(),
      mappings.getCoveredSpecies().end());
  for (const auto &species: sortedSpecies) {
    unsigned int spid = speciesToSpid.size();
    speciesToSpid.insert({species, spid});
  }
  for (auto geneTree: geneTrees) {
    for (auto leaf: geneTree->getLeaves()) {
      auto spid = speciesToSpid[leaf->label];
      leaf->data = (void *)(size_t)(spid);
    }
  }
}

using Split = BitVector;

static size_t getSpid(corax_unode_t *node)
{
  return (size_t)(node->data);
}

void getTaxonSplitsAux(corax_unode_t *node,
    Split &split)
{
  if (!node->next) {
    split.set(getSpid(node));
  } else {
    getTaxonSplitsAux(PLLUnrootedTree::getLeft(node), split);
    getTaxonSplitsAux(PLLUnrootedTree::getRight(node), split);
  }
  
}

corax_unode_t *getLCA(corax_unode_t *subtree,
    const BitVector &mask)
{
  if (!subtree->next) {
    auto res = mask.get(getSpid(subtree)) ? subtree : nullptr;
    return res;
  }
  auto lca1 = getLCA(PLLUnrootedTree::getLeft(subtree), mask);
  auto lca2 = getLCA(PLLUnrootedTree::getRight(subtree), mask);
  if (!lca1) {
    return lca2;
  } else if (!lca2) {
    return lca1;
  } else {
    return subtree;
  }
}


/**
 *  Locate the taxon with data field == spid
 *  Then traverses the right (resp. left)  subtree, fill
 *  a bitvector with 1s at the positions of the leaf data
 *  fields, and add it to leftSplits (resp. rightSplits)
 *
 *  Nothing happens if the taxon is not in the gene tree
 */
void getTaxonSplits(unsigned int spid,
    PLLUnrootedTree &geneTree,
    unsigned int speciesNumber,
    const BitVector &mask,
    std::vector<Split> &leftSplits,
    std::vector<Split> &rightSplits)
{
  // find the taxon
  corax_unode_t *taxon = nullptr;
  for (auto leaf: geneTree.getLeaves()) {
    if (getSpid(leaf) == spid) {
      taxon = leaf;
      break;
    }
  }
  if (!taxon) {
    return;
  }
  auto splitRoot = getLCA(taxon->back, mask);
  if (!splitRoot || !splitRoot->next) {
    // taxon is not present or only with one single taxon
    // from the current species tree
    return;
  }
  // fill the splits
  leftSplits.push_back(Split(speciesNumber, 0));
  rightSplits.push_back(Split(speciesNumber, 0));
  getTaxonSplitsAux(PLLUnrootedTree::getLeft(splitRoot), 
      leftSplits.back());
  getTaxonSplitsAux(PLLUnrootedTree::getRight(splitRoot), 
      rightSplits.back());
}


void getSplitAux(Node *node,
    StringToUint &speciesToSpid, 
    Split &split)
{
  if (!node->next) {
    split.set(speciesToSpid[node->label]);
  } else {
    getSplitAux(node->next->back, speciesToSpid, split);
    getSplitAux(node->next->next->back, speciesToSpid, split);
  } 
}

Split getSplit(Node *subtree,
    StringToUint &speciesToSpid, 
    unsigned int speciesNumber)
{
  Split split(speciesNumber, 0);
  getSplitAux(subtree, speciesToSpid, split);
  return split;
}

double interSize(const Split &sp1, const Split &sp2)
{
  auto inter = sp1;
  inter &=  sp2;
  auto count = inter.count();
  return double(count);
}

double computeDistance(const Split &b1,
    const Split &b2,
    const std::vector<Split> &leftSplits,
    const std::vector<Split> &rightSplits)
{
  double d = 0.0;
  for (unsigned int i = 0; i < leftSplits.size(); ++i) {
    const auto &s1 = leftSplits[i];
    const auto &s2 = rightSplits[i];
    auto v1 = interSize(b1, s1) + interSize(b2, s2); 
    auto v2 = interSize(b1, s2) + interSize(b2, s1); 
    auto diff = v1 * v2;
    d += diff * diff;
  }
  return d;
}

void insertTaxon(unsigned int spid,
    const std::string &label,
    StringToUint &speciesToSpid,
    unsigned int speciesNumber,
    const BitVector &mask,
    GeneTrees &perCoreGeneTrees,
    StepwiseTree &tree)
{
  std::vector<Split> leftSplits;
  std::vector<Split> rightSplits;
  for (auto geneTree: perCoreGeneTrees) {
    getTaxonSplits(spid, 
        *geneTree, 
        speciesNumber,
        mask,
        leftSplits,
        rightSplits);
  }
  Node *bestBranch = nullptr;
  double bestDistance = std::numeric_limits<double>::max();
  for (auto branch: tree.getBranches()) {
    auto leftBranchSplit = getSplit(branch, speciesToSpid, speciesNumber);
    auto rightBranchSplit = getSplit(branch->back, speciesToSpid, speciesNumber);
    auto d = computeDistance(leftBranchSplit,
        rightBranchSplit,
        leftSplits,
        rightSplits);
    if (d <= bestDistance) {
      bestDistance = d;
      bestBranch = branch;
    }
  }
  assert(bestBranch);
  //auto branch = tree.getRandomBranch();
  tree.addLeaf(label, bestBranch, spid);
}

void computeCoverages(GeneTrees &geneTrees,
    unsigned int speciesNumber,
    Coverages &coverages)
{
  for (auto &geneTree: geneTrees) {
    coverages.push_back(BitVector(speciesNumber, false));
    for (auto leaf: geneTree->getLeaves()) {
      coverages.back().set(getSpid(leaf));
    }
  }
}

void reorderTaxa(std::vector<std::string> &taxa,
    const Coverages &coverages,
    const StringToUint &speciesToSpid,
    bool random)
{
  if (random) {
    std::shuffle(taxa.begin(), taxa.end(), Random::getRNG());
    return;
  }
  unsigned int taxonNumber = taxa.size();
  std::unordered_set<std::string> remainingTaxa(taxa.begin(), taxa.end());
  auto firstTaxon = taxa[Random::getInt(taxa.size())]; 
  taxa.clear();
  BitVector mask(taxonNumber, false);
  
  taxa.push_back(firstTaxon);
  mask.set(speciesToSpid.at(firstTaxon));
  remainingTaxa.erase(firstTaxon);
  while (remainingTaxa.size()) {
    double bestScore = 0;
    std::string bestTaxon;
    for (auto taxon: remainingTaxa) {
      double score = 0;
      auto spid = speciesToSpid.at(taxon);
      double den = 0.0;
      for (auto &coverage: coverages) {
        if (coverage[spid]) {
          auto c = (coverage & mask).count();
          score += c * c;
          den += 1;
        }
      }
      //score = score / den;
      if (score > bestScore) {
        bestScore = score;
        bestTaxon = taxon;
      }
    }
    assert(bestTaxon.size());
    assert(speciesToSpid.find(bestTaxon) != speciesToSpid.end());
    Logger::info << "Adding " << bestTaxon << " with score " << bestScore << std::endl;
    auto bestSpid = speciesToSpid.at(bestTaxon);
    taxa.push_back(bestTaxon);
    mask.set(bestSpid);
    remainingTaxa.erase(bestTaxon);
  }
}





int main(int argc, char * argv[])
{
  // user-specified arguments
  Arguments arg(argc, argv);
  // all gene trees
  GeneTrees geneTrees;
  // trees assigned to the local parallel core/rank
  //GeneTrees perCoreGeneTrees;
  // map gene label to species label
  GeneSpeciesMapping mappings;
  // map species labels to species spid
  StringToUint speciesToSpid;
  // per family coverages
  Coverages coverages;


  init(arg);
  readGeneTrees(arg.inputGeneTreeFile, geneTrees);
  if (geneTrees.size() == 0) {
    Logger::error << "No input gene tree, aborting." << std::endl;
    exit(1);
  }
  extractMappings(arg, geneTrees, mappings, speciesToSpid); 
  unsigned int speciesNumber = mappings.getCoveredSpecies().size();
  auto speciesLabels = mappings.getCoveredSpecies();
  computeCoverages(geneTrees, speciesNumber, coverages);
  std::vector<std::string> taxonSequence(speciesLabels.begin(),
      speciesLabels.end());
  unsigned int iterations = 1;
  std::string outputFile = arg.prefix + ".newick";
  std::ofstream os(outputFile);
  for (unsigned int iter = 0; iter < iterations; ++iter) {
    reorderTaxa(taxonSequence, coverages, speciesToSpid, arg.randomInsertion);
    Logger::info << "RUnning search " << iter << std::endl; 
    BitVector mask(speciesNumber, false);
    StepwiseTree speciesTree(taxonSequence[0],
        taxonSequence[1],
        taxonSequence[2]);
    mask.set(speciesToSpid[taxonSequence[0]]);
    mask.set(speciesToSpid[taxonSequence[1]]);
    mask.set(speciesToSpid[taxonSequence[2]]);
    for (unsigned int taxonID = 3; taxonID < speciesNumber; ++taxonID) {
      auto label =  taxonSequence[taxonID];
      Logger::info << "inserting " << label << "..." << std::endl;
      auto spid = speciesToSpid[label];
      mask.set(spid);
      insertTaxon(spid,
          label,
          speciesToSpid, 
          speciesNumber, 
          mask,
          geneTrees,
          speciesTree);
     // Logger::info << speciesTree.getNewickString() << std::endl;
    }
    Logger::timed << "Saving tree into " << outputFile << std::endl;
    os << speciesTree.getNewickString() << std::endl;
  }
  return 0;

}
