#include <iostream>
#include <corax/corax.h>

#include "AsteroidOptimizer.hpp"
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <maths/Random.hpp>
#include "Arguments.hpp"
#include <trees/PLLUnrootedTree.hpp>
#include <trees/SplitHashtable.hpp>
#include <DistanceMethods/InternodeDistance.hpp>
#include <memory>
#include <set>
#include <algorithm>
#include <IO/ParallelOfstream.hpp>
#include <maths/bitvector.hpp>

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
 *  Creates a matrix of doubles of size NxM and fills it
 *  with value
 */
static DistanceMatrix initDistancetMatrix(unsigned int N,
    unsigned int M,
    double value)
{
  std::vector<double> nullDistances(M, value);
  return DistanceMatrix(N, nullDistances); 
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
  Logger::timed << "Number of gene trees:" << geneTrees.size() << std::endl;
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
    speciesToSpid.insert({species, speciesToSpid.size()});
  }
  Logger::timed << "Number of species: " 
    << speciesToSpid.size() 
    << std::endl;
}



/**
 *  Computes gidToSpid, which is a mapping between
 *  the gene IDs and the species IDs (for the gene 
 *  and species leaves). The gid is 
 *  also stored in the data field of the gene leaves
 */
void fillGidToSpid(std::vector<GeneCell> &geneCells,
    GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid)
{
  StringToUint labelToGid;
  for (auto &geneCell: geneCells) {
    auto &geneTree = *geneCell.geneTree;
    // the following line 
    // ensures that two gene trees with the same
    // leaf set also have the same leaf node indices:
    geneTree.sortLeaves();
    for (auto leaf: geneTree.getLeaves()) {
      intptr_t gid = leaf->node_index;
      const auto &species = mappings.getSpecies(leaf->label);
      auto spid = speciesToSpid.at(species);
      leaf->data = reinterpret_cast<void*>(gid);
      geneCell.gidToSpid.push_back(spid);
    }
  }
}

void computeGeneDistances(const Arguments &arg,
    std::vector<GeneCell *> &geneCells)
{
  Logger::timed << "Computing internode distances from the gene trees..." << std::endl;
  IDParam idParams;
  idParams.minBL = arg.minBL;
  for (auto cell: geneCells) {
    auto geneNumber = cell->geneTree->getLeavesNumber();
    cell->distanceMatrix = initDistancetMatrix(geneNumber,
        geneNumber,
        std::numeric_limits<double>::infinity());
    InternodeDistance::computeFromGeneTree(*cell->geneTree,
        cell->distanceMatrix, 
        idParams);
  }
}
  
/**
 *  Compute the coverage pattern of each gene tree
 */
void computeFamilyCoverage(
    std::vector<GeneCell> &geneCells,
    unsigned int speciesNumber)
{
  Logger::timed << "Computing coverage..." << std::endl;
  size_t total = 0;
  for (auto &geneCell: geneCells) {
    geneCell.coverage = BitVector(speciesNumber, false);
    for (auto geneLeaf: geneCell.geneTree->getLeaves()) {
      auto gid = reinterpret_cast<intptr_t>(geneLeaf->data);
      auto spid = geneCell.gidToSpid[gid];
      if (!geneCell.coverage[spid]) {
        total++;
        geneCell.coverage.set(spid);
      }
    }
  }
  Logger::timed << "Total number of genes: " << total << std::endl;
  Logger::timed << "Sampling proportion: " 
    << double(total) / double(geneCells.size() * speciesNumber)
    << std::endl;
} 

/**
 *  Assign to each MPI rank a subset of the gene 
 *  trees in order to parallelize computations
 */
void getPerCoreGeneCells(std::vector<GeneCell *> &geneCells,
    std::vector<GeneCell *> &perCoreGeneCells)
{
  auto K = geneCells.size();
  auto begin = ParallelContext::getBegin(K);
  auto end = ParallelContext::getEnd(K);
  perCoreGeneCells.clear();
  for (unsigned int k = begin; k < end; ++k) {
    perCoreGeneCells.push_back(geneCells[k]);
  }
}
    
void getDataNoCorrection(const std::vector<GeneCell *> &geneCells,
        std::vector<DistanceMatrix> &distanceMatrices,
        std::vector<BitVector> &coverages,
        UIntMatrix &gidToSpid,
        unsigned int speciesNumber)
{
  DistanceMatrix astridDistanceMatrix = initDistancetMatrix(speciesNumber,
      speciesNumber,
      0.0);
  DistanceMatrix denominator = initDistancetMatrix(speciesNumber,
      speciesNumber,
      0.0);
  for (auto cell: geneCells) {
    auto geneNumber = cell->gidToSpid.size();
    for (unsigned int igid = 0; igid < geneNumber; ++igid) {
      auto ispid = cell->gidToSpid[igid];
      for (unsigned int jgid = 0; jgid < geneNumber; ++jgid) {
        auto jspid = cell->gidToSpid[jgid];
        astridDistanceMatrix[ispid][jspid] += 
          cell->distanceMatrix[igid][jgid];
        denominator[ispid][jspid] += 1.0; 
      }
    }
  }
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      ParallelContext::sumDouble(astridDistanceMatrix[i][j]);
      ParallelContext::sumDouble(denominator[i][j]);
      if (denominator[i][j] != 0.0) {
        astridDistanceMatrix[i][j] /= denominator[i][j];
      }
    }
  }
  distanceMatrices.push_back(astridDistanceMatrix);
  coverages.push_back(BitVector (speciesNumber, true));
  gidToSpid.push_back(std::vector<unsigned int>());
  for (unsigned int i = 0;  i < speciesNumber; ++i) {
    gidToSpid.back().push_back(i);
  }
}

void getDataWithCorrection(const std::vector<GeneCell *> &geneCells,
        std::vector<DistanceMatrix> &distanceMatrices,
        std::vector<BitVector> &coverages,
        UIntMatrix &gidToSpid,
        unsigned int speciesNumber)
{
  DistanceMatrix astridDistanceMatrix = initDistancetMatrix(speciesNumber,
      speciesNumber,
      0.0);
  BitVector lastCoverage;
  for (auto cell: geneCells) {
    auto geneNumber = cell->gidToSpid.size();
    if (cell->coverage != lastCoverage) {
      distanceMatrices.push_back(cell->distanceMatrix);
      coverages.push_back(cell->coverage);
      gidToSpid.push_back(cell->gidToSpid);
    } else {
      auto &d = distanceMatrices.back();
      for (unsigned int igid = 0; igid < geneNumber; ++igid) {
        for (unsigned int jgid = 0; jgid < geneNumber; ++jgid) {
          d[igid][jgid] += cell->distanceMatrix[igid][jgid];
        }
      }
    }
    lastCoverage = cell->coverage;
  }
}

double optimize(PLLUnrootedTree &speciesTree,
    const std::vector<GeneCell *>  &geneCells,
    bool noCorrection)
{
  std::vector<GeneCell *> compressedCells;
  std::vector<DistanceMatrix> distanceMatrices;
  std::vector<BitVector> coverages;
  UIntMatrix gidToSpid;
  if (noCorrection) {
    getDataNoCorrection(geneCells,
        distanceMatrices,
        coverages,
        gidToSpid,
        speciesTree.getLeavesNumber());
  } else {
    getDataWithCorrection(geneCells,
        distanceMatrices,
        coverages,
        gidToSpid,
        speciesTree.getLeavesNumber());
  }
  Logger::timed << "Initializing optimizer with " << distanceMatrices.size() << " coverage patterns from " << geneCells.size() << " gene trees... " << std::endl;
  AsteroidOptimizer optimizer(speciesTree,
      coverages,
      gidToSpid,
      distanceMatrices);
  Logger::timed << "Starting tree search... " << std::endl;
  return optimizer.optimize();
}


ScoredTrees search(SpeciesTrees &startingSpeciesTrees,
  const std::vector<GeneCell *> &geneCells,
  bool noCorrection,
  int samples = -1)
{
  ScoredTrees scoredTrees;
  // Perform search
  for (auto speciesTree: startingSpeciesTrees) {
    ScoredTree st;
    st.tree = speciesTree;
    if (samples == -1) {
      st.score = optimize(*speciesTree, 
          geneCells, 
          noCorrection);
    } else {
      //assert(false); // todo re-enable!
      /*
      std::vector<DistanceMatrix> sampledDistanceMatrices;
      std::vector<BitVector> sampledPerFamilyCoverage;
      for (int i = 0; i < samples; ++i) {
        unsigned int s = Random::getInt(samples);
        sampledDistanceMatrices.push_back(asteroidDistanceMatrices[s]);
        sampledPerFamilyCoverage.push_back(asteroidPerFamilyCoverage[s]);
      }
      st.score = optimize(*speciesTree,
          sampledDistanceMatrices,
          sampledPerFamilyCoverage,
          noCorrection);
      */
    }
    scoredTrees.push_back(st);
  }
  std::sort(scoredTrees.begin(), scoredTrees.end());
  return scoredTrees;
}

TreePtr generateRandomSpeciesTree(std::unordered_set<std::string> &labels,
  const StringToUint &speciesToSpid)
{
  auto speciesTree = std::make_shared<PLLUnrootedTree>(
      labels);
  speciesTree->reindexLeaves(speciesToSpid);
  for (auto leaf: speciesTree->getLeaves()) {
    intptr_t spid = static_cast<intptr_t>(leaf->node_index);
    leaf->data = reinterpret_cast<void*>(spid);
  }
  return speciesTree;
}
  
Trees generateRandomSpeciesTrees(std::unordered_set<std::string> &labels,
    const StringToUint &speciesToSpid,
    unsigned int number)
{
  Trees trees;
  for (unsigned int i = 0; i < number; ++i) {
    auto speciesTree = generateRandomSpeciesTree(labels, 
        speciesToSpid);
    trees.push_back(speciesTree);
  } 
  return trees;
}



int main(int argc, char * argv[])
{
  // user-specified arguments
  Arguments arg(argc, argv);
  // all gene trees
  GeneTrees geneTrees;
  // trees assigned to the local parallel core/rank
  GeneTrees perCoreGeneTrees;
  // map gene label to species label
  GeneSpeciesMapping mappings;
  // map species labels to species spid
  StringToUint speciesToSpid;
  // per-family distance matrices
  std::vector<DistanceMatrix> distanceMatrices;
  // per family species coverage
  std::vector<BitVector> perFamilyCoverage;
   
  UIntMatrix gidToSpid;

  init(arg);
  readGeneTrees(arg.inputGeneTreeFile, geneTrees);
  extractMappings(arg, geneTrees, mappings, speciesToSpid); 
  unsigned int speciesNumber = mappings.getCoveredSpecies().size();
  auto speciesLabels = mappings.getCoveredSpecies();
 

  std::vector<GeneCell> allGeneCells;
  for (auto geneTree: geneTrees) {
    allGeneCells.push_back(GeneCell(geneTree));
  }
  fillGidToSpid(allGeneCells, mappings, speciesToSpid);
  computeFamilyCoverage(allGeneCells, speciesNumber);
  std::unordered_map<BitVector, std::vector<GeneCell *> > cellMap;
  for (auto &cell: allGeneCells) {
    if (cellMap.end() == cellMap.find(cell.coverage)) {
      cellMap.insert({cell.coverage, std::vector<GeneCell *>()});
    }
    cellMap[cell.coverage].push_back(&cell);
  }
  std::vector<GeneCell *> geneCells;
  std::vector<GeneCell *> perCoreGeneCells;
  for (auto it: cellMap) {
    for (auto cell: it.second) {
      geneCells.push_back(cell);
    }
  }
  getPerCoreGeneCells(geneCells, perCoreGeneCells);
  computeGeneDistances(arg, perCoreGeneCells);

  SpeciesTrees startingSpeciesTrees =
    generateRandomSpeciesTrees(speciesLabels, 
        speciesToSpid,
        std::max(arg.randomStartingTrees, 1u));
  if (arg.randomStartingTrees < 1) {
    optimize(*startingSpeciesTrees[0], 
          perCoreGeneCells,
          true);
  }
  ScoredTrees scoredTrees = search(startingSpeciesTrees,
      perCoreGeneCells,
      arg.noCorrection);
  std::string outputSpeciesTreeFile = arg.prefix + ".bestTree.newick";
  std::string outputSpeciesTreesFile = arg.prefix + ".allTrees.newick";
  std::string outputScoresFile = arg.prefix + ".scores.txt"; 
  Logger::timed << "Final tree scores: " << std::endl;
  Logger::timed << "Printing the best tree into " << outputSpeciesTreeFile << std::endl;
  auto & bestTree = *scoredTrees[0].tree;
  bestTree.save(outputSpeciesTreeFile);
  ParallelOfstream treesOs(outputSpeciesTreesFile);
  ParallelOfstream scoresOs(outputScoresFile);
  for (const auto &st: scoredTrees) {
    treesOs << st.tree->getNewickString() << std::endl;
    scoresOs << st.score << std::endl;
    Logger::info << st.score << std::endl;
  }


  if (arg.bootstrapTrees > 0) {
    assert(false);
    /*
    SpeciesTrees startingBSTrees =
      generateRandomSpeciesTrees(speciesLabels, 
          speciesToSpid,
          arg.bootstrapTrees);
    ScoredTrees finalBSTrees = search(startingBSTrees,
        distanceMatrices,
        perFamilyCoverage,
        arg.noCorrection,
        perCoreGeneTrees.size());
    std::string bsTreesFile = arg.prefix + ".bsTrees.newick";
    ParallelOfstream bsOs(bsTreesFile);
    SplitHashtable splits;
    for (const auto &st: finalBSTrees) {
      bsOs << st.tree->getNewickString() << std::endl;
      splits.addTree(*st.tree);
    }
    splits.computeSupportValues(bestTree);
    bestTree.save(outputSpeciesTreeFile);
    */
  } 
  

  close();
  return 0;

}
