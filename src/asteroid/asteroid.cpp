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

void setGeneTreeSpid(PLLUnrootedTree &geneTree,
    const GeneSpeciesMapping &mappings,
    const StringToUint &speciesToSpid)
{
  for (auto leaf: geneTree.getLeaves()) {
    const auto &species = mappings.getSpecies(leaf->label);
    intptr_t spid = static_cast<intptr_t>(speciesToSpid.at(species));
    leaf->data = reinterpret_cast<void*>(spid);
  }
}

static DistanceMatrix initDistancetMatrix(unsigned int N,
    unsigned int M,
    double value)
{
  std::vector<double> nullDistances(M, value);
  return DistanceMatrix(N, nullDistances); 
}
  
void fillDistanceMatricesAstrid(GeneTrees &geneTrees, 
    unsigned int speciesNumber, 
    const IDParam &params,
    std::vector<DistanceMatrix> &distanceMatrices)
{
  DistanceMatrix sum = initDistancetMatrix(speciesNumber,
      speciesNumber, 
      0.0);
  DistanceMatrix denominator = initDistancetMatrix(speciesNumber,
      speciesNumber, 
      0.0);
  for (unsigned int k = 0; k < geneTrees.size(); ++k) {
  DistanceMatrix temp = initDistancetMatrix(speciesNumber,
      speciesNumber, 
      std::numeric_limits<double>::infinity());
    auto &geneTree = *geneTrees[k];
    InternodeDistance::computeFromGeneTree(geneTree, temp, params);
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      for (unsigned int j = 0; j < speciesNumber; ++j) {
        if (temp[i][j] != std::numeric_limits<double>::infinity()) {
          sum[i][j] += temp[i][j];
          temp[i][j] = std::numeric_limits<double>::infinity();
          denominator[i][j] += 1.0;
        }
      }
    }
  }
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      if (i ==j) {
        continue;
      }
      if (denominator[i][j] == 0.0) {
        Logger::info << "Warning: missing entry in the distance matrix" << std::endl;
        
      } else {
        sum[i][j] /= denominator[i][j];
      }
    }
  }
  distanceMatrices.push_back(sum);
}

void fillDistanceMatricesAsteroid(GeneTrees &geneTrees, 
    unsigned int speciesNumber, 
    const IDParam &params,
    std::vector<DistanceMatrix> &distanceMatrices)
{
  distanceMatrices.resize(geneTrees.size());
  for (unsigned int k = 0; k < geneTrees.size(); ++k) {
    auto &d = distanceMatrices[k];
    auto &geneTree = *geneTrees[k];
    d = initDistancetMatrix(speciesNumber,
        speciesNumber,
        std::numeric_limits<double>::infinity());
    InternodeDistance::computeFromGeneTree(geneTree, d, params);
  }
}

void init(Arguments &arg)
{
  ParallelContext::init(0);
  Logger::init();
  Logger::timed << "Starting Asteroid..." << std::endl;
  arg.printCommand();
  Random::setSeed(static_cast<unsigned int>(arg.seed));
}

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

void extractMappings(const Arguments &arg,
    GeneTrees &geneTrees,
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


void applyMappings(GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid,
    GeneTrees &geneTrees)
{
  for (auto geneTree: geneTrees) {
    setGeneTreeSpid(*geneTree, mappings, speciesToSpid);
  }
}

void computeGeneDistances(const Arguments &arg,
    GeneTrees &geneTrees,
    unsigned int speciesNumber,
    std::vector<DistanceMatrix> &distanceMatrices)
{
  Logger::timed << "Computing internode distances from the gene trees..." << std::endl;
  IDParam idParams;
  idParams.minBL = arg.minBL;
 
  fillDistanceMatricesAsteroid(geneTrees, 
      speciesNumber, 
      idParams, 
      distanceMatrices);
}
  
void computeFamilyCoverage(
    GeneTrees &geneTrees,
    unsigned int speciesNumber,
    BoolMatrix &perFamilyCoverage)
{
  Logger::timed << "Computing coverage..." << std::endl;
  perFamilyCoverage = BoolMatrix(geneTrees.size(),
      std::vector<bool>(speciesNumber, false));
  for (unsigned int k = 0; k < geneTrees.size(); ++k) {
    for (auto geneLeaf: geneTrees[k]->getLeaves()) {
      auto spid = reinterpret_cast<intptr_t>(geneLeaf->data);
      perFamilyCoverage[k][spid] = true;
    }
  }
} 

void getPerCoreGeneTrees(GeneTrees &geneTrees,
    GeneTrees &perCoreGeneTrees)
{
  auto K = geneTrees.size();
  auto begin = ParallelContext::getBegin(K);
  auto end = ParallelContext::getEnd(K);
  perCoreGeneTrees.clear();
  for (unsigned int k = begin; k < end; ++k) {
    perCoreGeneTrees.push_back(geneTrees[k]);
  }
}
    
void getDataNoCorrection(
        const std::vector<DistanceMatrix> &distanceMatrices,
        std::vector<DistanceMatrix> &astridDistanceMatrices,
        BoolMatrix & astridPerFamilyCoverage)
{
  unsigned int speciesNumber = distanceMatrices[0].size();
  unsigned int K = distanceMatrices.size(); 
  DistanceMatrix astridDistanceMatrix = initDistancetMatrix(speciesNumber,
      speciesNumber,
      0.0);
  DistanceMatrix denominator = initDistancetMatrix(speciesNumber,
      speciesNumber,
      0.0);
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      for (unsigned int k = 0; k < K; ++k) {
        if (std::numeric_limits<double>::infinity() 
            != distanceMatrices[k][i][j]) {
          astridDistanceMatrix[i][j] += distanceMatrices[k][i][j];
          denominator[i][j] += 1.0; 
        }
      }
      ParallelContext::sumDouble(astridDistanceMatrix[i][j]);
      ParallelContext::sumDouble(denominator[i][j]);
      if (denominator[i][j] != 0.0) {
        astridDistanceMatrix[i][j] /= denominator[i][j];
      }
    }
  }
  astridDistanceMatrices.push_back(astridDistanceMatrix);
  std::vector<bool> astridCoverage(speciesNumber, true);
  astridPerFamilyCoverage.push_back(astridCoverage);
}

double optimize(PLLUnrootedTree &speciesTree,
    const std::vector<DistanceMatrix>  &asteroidDistanceMatrices,
    const BoolMatrix &asteroidPerFamilyCoverage,
    bool noCorrection)
{
  const std::vector<DistanceMatrix> *distanceMatrices = nullptr;
  const BoolMatrix *perFamilyCoverage = nullptr;
  std::vector<DistanceMatrix> astridDistanceMatrices;
  BoolMatrix astridPerFamilyCoverage;
  if (noCorrection) {
    getDataNoCorrection(asteroidDistanceMatrices,
        astridDistanceMatrices,
        astridPerFamilyCoverage);
    perFamilyCoverage = &astridPerFamilyCoverage;
    distanceMatrices = &astridDistanceMatrices;
  } else {
    perFamilyCoverage = &asteroidPerFamilyCoverage;
    distanceMatrices = &asteroidDistanceMatrices;
  }
  Logger::timed << "Initializing optimizer... " << std::endl;
  AsteroidOptimizer optimizer(speciesTree,
      *perFamilyCoverage,
      *distanceMatrices);
  Logger::timed << "Starting tree search... " << std::endl;
  return optimizer.optimize();
}

ScoredTrees search(SpeciesTrees &startingSpeciesTrees,
  const std::vector<DistanceMatrix> &asteroidDistanceMatrices,
  const BoolMatrix &asteroidPerFamilyCoverage,
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
          asteroidDistanceMatrices, 
          asteroidPerFamilyCoverage,
          noCorrection);
    } else {
      std::vector<DistanceMatrix> sampledDistanceMatrices;
      BoolMatrix sampledPerFamilyCoverage;
      for (int i = 0; i < samples; ++i) {
        unsigned int s = Random::getInt(samples);
        sampledDistanceMatrices.push_back(asteroidDistanceMatrices[s]);
        sampledPerFamilyCoverage.push_back(asteroidPerFamilyCoverage[s]);
      }
      st.score = optimize(*speciesTree,
          sampledDistanceMatrices,
          sampledPerFamilyCoverage,
          noCorrection);

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
  Arguments arg(argc, argv);
  GeneTrees geneTrees;
  GeneTrees perCoreGeneTrees;
  StringToUint speciesToSpid;
  std::vector<DistanceMatrix> distanceMatrices;
  BoolMatrix perFamilyCoverage;
  
  init(arg);
  readGeneTrees(arg.inputGeneTreeFile, geneTrees);
  GeneSpeciesMapping mappings;
  extractMappings(arg, geneTrees, mappings, speciesToSpid); 
  unsigned int speciesNumber = mappings.getCoveredSpecies().size();
  auto speciesLabels = mappings.getCoveredSpecies();
 

  getPerCoreGeneTrees(geneTrees, perCoreGeneTrees);
  applyMappings(mappings, speciesToSpid, perCoreGeneTrees); 
  computeFamilyCoverage(perCoreGeneTrees, speciesNumber, perFamilyCoverage);
  computeGeneDistances(arg, perCoreGeneTrees, speciesNumber, distanceMatrices); 
  
   

  SpeciesTrees startingSpeciesTrees =
    generateRandomSpeciesTrees(speciesLabels, 
        speciesToSpid,
        arg.randomStartingTrees);
  ScoredTrees scoredTrees = search(startingSpeciesTrees,
      distanceMatrices,
      perFamilyCoverage,
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
  } 
  

  close();
  return 0;

}
