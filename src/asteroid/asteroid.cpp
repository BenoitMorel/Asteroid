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
#include <DistanceMethods/StepwiseAsteroid.hpp>
#include <memory>
#include <set>
#include <algorithm>
#include <IO/ParallelOfstream.hpp>
#include <maths/bitvector.hpp>

#include <DistanceMethods/AsteroidTypes.hpp>


/**
 *  Read all trees from inputGeneTreeFile
 */
static void parseTrees(const std::string &inputGeneTreeFile,
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
    if (geneTrees.size() % 100000 == 0) {
      Logger::info << geneTrees.size() << std::endl;
    }
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
static void init(Arguments &arg)
{
  
  ParallelContext::init(nullptr);
  Logger::init();
  Logger::timed << "Asteroid v1.0" << std::endl;
  arg.printCommand();
  Random::setSeed(static_cast<unsigned int>(arg.seed));
  assert(ParallelContext::isRandConsistent());
}

/**
 *  Closes the program state
 */
static void close()
{
  Logger::close();
  ParallelContext::finalize();
}

static void readGeneTrees(const std::string &inputGeneTreeFile,
    GeneTrees &geneTrees)
{
  Logger::timed << "Parsing gene trees..." << std::endl;
  parseTrees(inputGeneTreeFile,
      geneTrees);
  Logger::timed << "End of parsing gene trees" << std::endl;
}
  
static void readUserWeights(const std::string &userWeightsFile, 
    std::vector<double> &userWeights)
{
  Logger::timed << "Parsing user-specified weights..." << std::endl;
  std::ifstream reader(userWeightsFile);
  std::string weightStr;
  while (std::getline(reader, weightStr)) {
    if (weightStr.size() < 2) {
      continue;
    }
    double weight = atof(weightStr.c_str());
    userWeights.push_back(weight);
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
static void extractMappings(const Arguments &arg,
    const GeneTrees &geneTrees,
    GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid)
{
  Logger::timed << "Mapping gene to species..." << std::endl;
  if (!arg.geneSpeciesMapping.size()) {
    Logger::timed << "Extracting mappings from the gene trees..." << std::endl;
    for (auto geneTree: geneTrees) {
      mappings.fillFromGeneTree(*geneTree, true);
    }
  } else {
    Logger::timed << "Extracting mappings from " << arg.geneSpeciesMapping << std::endl;
    mappings.fillFromMappingFile(arg.geneSpeciesMapping); 
  }
  std::set<std::string> sortedSpecies(mappings.getCoveredSpecies().begin(),
      mappings.getCoveredSpecies().end());
  for (const auto &species: sortedSpecies) {
    speciesToSpid.insert({species, speciesToSpid.size()});
  }
}



/**
 *  Computes gidToSpid, which is a mapping between
 *  the gene IDs and the species IDs (for the gene 
 *  and species leaves). The gid is 
 *  also stored in the data field of the gene leaves
 */
static void fillGidToSpid(std::vector<GeneCell> &geneCells,
    GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid)
{
  StringToUint labelToGid;
  for (auto &geneCell: geneCells) {
    std::unordered_map<unsigned int, unsigned int> spidToGid;
    auto &geneTree = *geneCell.geneTree;
    // the following line 
    // ensures that two gene trees with the same
    // leaf set also have the same leaf node indices:
    geneTree.sortLeaves();
    for (auto leaf: geneTree.getLeaves()) {
      const auto &species = mappings.getSpecies(leaf->label);
      auto spid = speciesToSpid.at(species);
      intptr_t gid = 0;
      auto it = spidToGid.find(spid);
      if (it != spidToGid.end()) {
        gid = it->second;
      } else {
        gid = spidToGid.size();
        spidToGid.insert({spid, gid});
      }
      //intptr_t gid = leaf->node_index;
      leaf->data = reinterpret_cast<void*>(gid);
      //geneCell.gidToSpid.push_back(spid);
    }
    geneCell.gidToSpid.resize(spidToGid.size());
    for (auto it: spidToGid) {
      geneCell.gidToSpid[it.second] = it.first;
    }
  }
}

static void computeGeneDistances(const Arguments &arg,
    const std::vector<GeneCell *> &geneCells)
{
  Logger::timed << "Computing internode distances from the gene trees..." << std::endl;
  IDParam idParams;
  idParams.minBL = arg.minBL;
  idParams.useBL = arg.useGeneBL;
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
static void computeFamilyCoverage(
    std::vector<GeneCell> &geneCells,
    unsigned int speciesNumber)
{
  Logger::timed << "Computing coverage..." << std::endl;
  size_t total = 0;
  std::unordered_set<BitVector> coverages;
  for (auto &geneCell: geneCells) {
    geneCell.coverage = BitVector(speciesNumber, false);
    for (auto geneLeaf: geneCell.geneTree->getLeaves()) {
      auto gid = static_cast<unsigned int>(reinterpret_cast<intptr_t>(geneLeaf->data));
      auto spid = geneCell.gidToSpid[gid];
      if (!geneCell.coverage[spid]) {
        total++;
        geneCell.coverage.set(spid);
      }
    }
    coverages.insert(geneCell.coverage);
  }
  Logger::info << "\tNumber of species: " << speciesNumber << std::endl;
  Logger::info << "\tTotal number of gene trees: " << geneCells.size() << std::endl;
  Logger::info << "\tNumber of distinct coverage patterns:" << coverages.size() << std::endl;
  Logger::info << "\tTotal number of gene leaves: " << total << std::endl;
  Logger::info << "\tSampling proportion: " 
    << double(total) / double(geneCells.size() * speciesNumber)
    << std::endl;
}

static std::vector<unsigned int> getPerCoreUniformSampling(unsigned int K)
{
  auto begin = ParallelContext::getBegin(K);
  auto end = ParallelContext::getEnd(K);
  return std::vector<unsigned int>(end - begin, 1);
}

static std::vector<unsigned int> getPerCoreBSSampling(unsigned int K) 
{
  std::vector<unsigned int> sampling(K, 0);
  for (unsigned int i = 0; i < K; ++i) {
    sampling[Random::getUInt(K)]++;
  }
  std::vector<unsigned int> perCoreSampling;
  auto begin = ParallelContext::getBegin(K);
  auto end = ParallelContext::getEnd(K);
  for (unsigned int i = begin; i < end; ++i) {
    perCoreSampling.push_back(sampling[i]);
  }
  return perCoreSampling;
}


/**
 *  Assign to each MPI rank a subset of the gene 
 *  trees in order to parallelize computations
 */
static void getPerCoreGeneCells(std::vector<GeneCell *> &geneCells,
    std::vector<GeneCell *> &perCoreGeneCells)
{
  auto K = static_cast<unsigned int>(geneCells.size());
  auto begin = ParallelContext::getBegin(K);
  auto end = ParallelContext::getEnd(K);
  perCoreGeneCells.clear();
  for (unsigned int k = begin; k < end; ++k) {
    perCoreGeneCells.push_back(geneCells[k]);
  }
}
    
static void getDataNoCorrection(const std::vector<GeneCell *> &geneCells,
        const std::vector<unsigned int> &weights,
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
  for (unsigned int k = 0; k < geneCells.size(); ++k ) {
    if (!weights[k]) {
      continue;
    } 
    auto cell = geneCells[k];
    double weight = static_cast<double>(weights[k]);
    auto geneNumber = cell->gidToSpid.size();
    for (unsigned int igid = 0; igid < geneNumber; ++igid) {
      auto ispid = cell->gidToSpid[igid];
      for (unsigned int jgid = 0; jgid < geneNumber; ++jgid) {
        auto jspid = cell->gidToSpid[jgid];
        astridDistanceMatrix[ispid][jspid] += 
          cell->distanceMatrix[igid][jgid] * weight * cell->userWeight;
        denominator[ispid][jspid] += weight; 
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

static void getDataWithCorrection(const std::vector<GeneCell *> &geneCells,
        const std::vector<unsigned int> &weights,
        std::vector<DistanceMatrix> &distanceMatrices,
        std::vector<BitVector> &coverages,
        UIntMatrix &gidToSpid)
{
  BitVector lastCoverage;
  for (unsigned int k = 0; k < geneCells.size(); ++k) {
    if (!weights[k]) {
      continue;
    } 
    auto cell = geneCells[k];
    double weight = static_cast<double>(weights[k]);
    auto geneNumber = cell->gidToSpid.size();
    if (cell->coverage != lastCoverage) {
      distanceMatrices.push_back(initDistancetMatrix(geneNumber, geneNumber, 0.0));
      coverages.push_back(cell->coverage);
      gidToSpid.push_back(cell->gidToSpid);
    }
    auto &d = distanceMatrices.back();
    for (unsigned int igid = 0; igid < geneNumber; ++igid) {
      for (unsigned int jgid = 0; jgid < geneNumber; ++jgid) {
        d[igid][jgid] += cell->distanceMatrix[igid][jgid] * weight * cell->userWeight;
      }
    }
    lastCoverage = cell->coverage;
  }
}

static double optimize(PLLUnrootedTree &speciesTree,
    const std::vector<GeneCell *>  &geneCells,
    const std::vector<unsigned int> &weights,
    bool noCorrection,
    bool verbose)
{
  std::vector<DistanceMatrix> distanceMatrices;
  std::vector<BitVector> coverages;
  UIntMatrix gidToSpid;
  if (noCorrection) {
    getDataNoCorrection(geneCells,
        weights,
        distanceMatrices,
        coverages,
        gidToSpid,
        speciesTree.getLeavesNumber());
  } else {
    getDataWithCorrection(geneCells,
        weights,
        distanceMatrices,
        coverages,
        gidToSpid);
  }
  AsteroidOptimizer optimizer(speciesTree,
      coverages,
      gidToSpid,
      distanceMatrices,
      verbose);
  return optimizer.optimize();
}


static ScoredTrees search(const SpeciesTrees &startingSpeciesTrees,
  const std::vector<GeneCell *> &geneCells,
  const std::vector<unsigned int> &weights,
  bool noCorrection,
  bool verbose,
  unsigned int &index) 
{
  ScoredTrees scoredTrees;
  // Perform search
  for (auto speciesTree: startingSpeciesTrees) {
    Logger::timed << "Starting tree search " << index++ << std::endl;
    ScoredTree st;
    st.tree = speciesTree;
    st.score = optimize(*speciesTree, 
          geneCells, 
          weights,
          noCorrection,
          verbose);
    Logger::info << "  score: " << st.score << std::endl;
    scoredTrees.push_back(st);
  }
  std::sort(scoredTrees.begin(), scoredTrees.end());
  return scoredTrees;
}

static TreePtr generateRandomSpeciesTree(std::unordered_set<std::string> &labels,
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
  
static Trees generateRandomSpeciesTrees(std::unordered_set<std::string> &labels,
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

static void reorderTaxa(std::vector<std::string> &taxa,
    const std::vector<GeneCell *> &geneCells,
    const StringToUint &speciesToSpid,
    bool random)
{
  if (random) {
    std::shuffle(taxa.begin(), taxa.end(), Random::getRNG());
    return;
  }
  auto taxonNumber = taxa.size();
  std::unordered_set<std::string> remainingTaxa(taxa.begin(), taxa.end());
  auto firstTaxon = taxa[Random::getUInt(static_cast<unsigned int>(taxa.size()))]; 
  taxa.clear();
  BitVector mask(taxonNumber, false);
 
  Logger::info << "Adding " << firstTaxon << std::endl;
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
      for (auto cell: geneCells) {
        const auto &coverage = cell->coverage;
        if (coverage[spid]) {
          auto c = (coverage & mask).count();
          score += c * c;
          den += 1;
        }
      }
      ParallelContext::sumDouble(score);
      ParallelContext::sumDouble(den);
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

static void generateStepwiseTree(const StringToUint &speciesToSpid,
    const std::vector<GeneCell *> &geneCells,
    const std::string &outputTree)
{
  std::vector<std::string> labels;
  for (const auto &pair: speciesToSpid) {
    labels.push_back(pair.first);
  }
  reorderTaxa(labels, geneCells, speciesToSpid, false);
  StepwiseAsteroid stepwise(speciesToSpid,
      geneCells);
  for (const auto &label: labels) {
    stepwise.insertLabel(label);
  }
  stepwise.exportTree(outputTree);
}


struct AsteroidInstance {
  GeneTrees geneTrees;
  std::vector<GeneCell> allGeneCells;
  std::vector<GeneCell *> geneCells;
  std::vector<GeneCell *> perCoreGeneCells;
  StringToUint speciesToSpid;
  // trees assigned to the local parallel core/rank
  GeneTrees perCoreGeneTrees;
  std::unordered_set<std::string> speciesLabels;
};


void initInstance(const Arguments &arg,
    const std::string &inputGeneTreeFile,
    AsteroidInstance &instance)
{
  // user-given gene tree weights
  std::vector<double> userWeights;
  // map gene label to species label
  GeneSpeciesMapping mappings;
  
  readGeneTrees(inputGeneTreeFile, instance.geneTrees);
  if (instance.geneTrees.size() == 0) {
    Logger::info << "Please provide at least one input gene tree" << std::endl;
    ParallelContext::abort(1);
  }
  extractMappings(arg, 
      instance.geneTrees, 
      mappings, 
      instance.speciesToSpid); 
  auto speciesNumber = static_cast<unsigned int>(mappings.getCoveredSpecies().size());
  instance.speciesLabels = mappings.getCoveredSpecies();
  for (auto geneTree: instance.geneTrees) {
    instance.allGeneCells.push_back(GeneCell(geneTree));
  }

  if (arg.inputWeights.size()) {
    assert(false);
    readUserWeights(arg.inputWeights, userWeights);
    if (userWeights.size() != instance.allGeneCells.size()) {
      Logger::info << "Error: the number of weights (" << userWeights.size() << 
        ") is not equal to the number of input gene trees (" << instance.allGeneCells.size() <<
        "). Aborting. " << std::endl;
    }
    for (unsigned int k = 0; k < instance.allGeneCells.size(); ++k) {
      instance.allGeneCells[k].userWeight = userWeights[k];
    }
  }
  fillGidToSpid(instance.allGeneCells, mappings, instance.speciesToSpid);
  computeFamilyCoverage(instance.allGeneCells, speciesNumber);
  std::unordered_map<BitVector, std::vector<GeneCell *> > cellMap;

  
  for (auto &cell: instance.allGeneCells) {
    if (cell.coverage.count() < 4) {
      continue;
    }
    if (cellMap.end() == cellMap.find(cell.coverage)) {
      cellMap.insert({cell.coverage, std::vector<GeneCell *>()});
    }
    cellMap[cell.coverage].push_back(&cell);
  }
  
  for (auto it: cellMap) {
    instance.geneCells.insert(instance.geneCells.end(), 
        it.second.begin(), 
        it.second.end());
  }
  getPerCoreGeneCells(instance.geneCells, instance.perCoreGeneCells);
  computeGeneDistances(arg, instance.perCoreGeneCells);
}


int main(int argc, char * argv[])
{
  // user-specified arguments
  Arguments arg(argc, argv);
  init(arg);
  
  if (!arg.isOk()) {
    Logger::info << "Aborting" << std::endl;
    close();
    return 1;
  }


  AsteroidInstance instance;
  initInstance(arg,
      arg.inputGeneTreeFile,
      instance);



  std::string outputSpeciesTreeFile = arg.prefix + ".bestTree.newick";
  if (arg.stepwise) {
    generateStepwiseTree(instance.speciesToSpid, 
        instance.perCoreGeneCells,
        outputSpeciesTreeFile);
    close();
    return 0;
  }

  auto K = static_cast<unsigned int>(instance.geneCells.size());
  SpeciesTrees startingSpeciesTrees =
    generateRandomSpeciesTrees(instance.speciesLabels, 
        instance.speciesToSpid,
        std::max(arg.randomStartingTrees, 1u));
  if (arg.randomStartingTrees < 1) {
    optimize(*startingSpeciesTrees[0], 
          instance.perCoreGeneCells,
          getPerCoreUniformSampling(K),
          true,
          false);
  }
  Logger::info << std::endl;
  unsigned int index = 0;
  Logger::timed << "Asteroid will now search for the best scoring tree from " << startingSpeciesTrees.size() << " starting trees" << std::endl;
  ScoredTrees scoredTrees = search(startingSpeciesTrees,
      instance.perCoreGeneCells,
      getPerCoreUniformSampling(K),
      arg.noCorrection,
      startingSpeciesTrees.size() == 1,
      index);
  std::string outputSpeciesTreesFile = arg.prefix + ".allTrees.newick";
  std::string outputScoresFile = arg.prefix + ".scores.txt"; 
  Logger::timed << "Printing the best tree into " << outputSpeciesTreeFile << std::endl;
  auto & bestTree = *scoredTrees[0].tree;
  bestTree.setAllBranchLengths(1.0);
  bestTree.save(outputSpeciesTreeFile);
  ParallelOfstream treesOs(outputSpeciesTreesFile);
  ParallelOfstream scoresOs(outputScoresFile);
  for (const auto &st: scoredTrees) {
    treesOs << st.tree->getNewickString() << std::endl;
    scoresOs << st.score << std::endl;
  }
  if (arg.bootstrapReplicates > 0) {
    Logger::info << std::endl;
    Logger::timed << "Asteroid will now compute " << arg.bootstrapReplicates << " bootstrap trees" << std::endl;
    assert(arg.inputBSGeneTreeFile.size());
    AsteroidInstance bsInstance;
    initInstance(arg,
        arg.inputBSGeneTreeFile,
        bsInstance);

    std::string bsTreesFile = arg.prefix + ".bsTrees.newick";
    ParallelOfstream bsOs(bsTreesFile);
    SplitHashtable splits;
    index = 0;
    for (unsigned int i = 0; i < arg.bootstrapReplicates; ++i) {
      SpeciesTrees startingBSTrees =
        generateRandomSpeciesTrees(bsInstance.speciesLabels, 
            bsInstance.speciesToSpid,
            1);
      ScoredTrees scoredBSTree = search(startingBSTrees,
        bsInstance.perCoreGeneCells,
        getPerCoreBSSampling(K),
        arg.noCorrection,
        false,
        index);
      bsOs << scoredBSTree[0].tree->getNewickString() << std::endl;
      splits.addTree(*scoredBSTree[0].tree);
    }
    splits.computeSupportValues(bestTree);
    bestTree.setAllBranchLengths(1.0);
    bestTree.save(outputSpeciesTreeFile);
    Logger::timed << "Updating the best-scoring tree with support values into " << outputSpeciesTreeFile << std::endl;
  } 
  close();
  return 0;

}
