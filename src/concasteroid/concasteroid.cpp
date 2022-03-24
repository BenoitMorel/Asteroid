#include <iostream>
#include <corax/corax.h>

#include <DistanceMethods/AsteroidOptimizer.hpp>
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

void parseDistanceMatrices(const std::string &distanceMatrixFile,
    GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid,
    std::vector<BitVector> &perFamilyCoverage,
    UIntMatrix &perFamilyGidToSpid,
    std::vector<DistanceMatrix> &distanceMatrices)
{
  auto speciesNumber = mappings.getCoveredSpecies().size();
  // Read all trees
  std::ifstream fileListReader(distanceMatrixFile);
  std::string line;
  unsigned int validCount = 0;
  unsigned int lines = 0;
  while (std::getline(fileListReader, line)) {
    lines++;
    std::ifstream is(line);
    if (is.peek() == std::ifstream::traits_type::eof()) {
      Logger::error << "Wrning: empty distance matrix" << std::endl;
      continue;
    }
    unsigned int N; // matrix size
    is >> N;
    auto coverage = BitVector(speciesNumber, false);
    auto matrix = DistanceMatrix(N, std::vector<double>(N));
    auto gidToSpid = std::vector<unsigned int>(N);
    bool valid = true;
    for (unsigned int i = 0; i < N && valid; ++i) {
      std::string gene;
      is >> gene;
      auto species = mappings.getSpecies(gene);
      auto spid = speciesToSpid[species];
      coverage.set(spid);
      gidToSpid[i] = spid;
      for (unsigned int j = 0; j < N; ++j) {
        std::string temp;
        is >> temp;
        if (temp == "NA") {
          valid = false;
          break;
        }
        matrix[i][j] = std::stod(temp);
        /*
        if (matrix[i][j] > 4.0) {
          valid = false;
          break;
        }
        */
      }
    }
    if (valid) {
      perFamilyCoverage.push_back(coverage);
      perFamilyGidToSpid.push_back(gidToSpid);
      distanceMatrices.push_back(matrix);
      validCount++;
    } else {
      Logger::info << "Warning: invalid distance matrix" << std::endl;
    }
  }
  Logger::timed << "Valid entries: " << validCount << "/" << lines << std::endl; 
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

/**
 *  Fills the gene-to-species mapping object:
 *  - from the input mapping file
 *  Identifies all the species labels and assign
 *  them a SPID
 *
 */
void extractMappings(const Arguments &arg,
    GeneSpeciesMapping &mappings,
    StringToUint &speciesToSpid)
{
  Logger::timed << "Mapping gene to species..." << std::endl;
    mappings.fillFromMappingFile(arg.geneSpeciesMapping); 
  std::set<std::string> sortedSpecies(mappings.getCoveredSpecies().begin(),
      mappings.getCoveredSpecies().end());
  for (const auto &species: sortedSpecies) {
    speciesToSpid.insert({species, speciesToSpid.size()});
  }
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
  // map gene label to species label
  GeneSpeciesMapping mappings;
  // map species labels to species spid
  StringToUint speciesToSpid;
  std::vector<std::string> spidToSpecies;
  // per-family distance matrices
  std::vector<DistanceMatrix> distanceMatrices;
  // per family species coverage
  std::vector<BitVector> perFamilyCoverage;
  // per family gid to spid
  UIntMatrix perFamilyGidToSpid;
   

  init(arg);
  extractMappings(arg, mappings, speciesToSpid); 
  
  parseDistanceMatrices(arg.distanceMatrixListFile, 
      mappings,
      speciesToSpid,
      perFamilyCoverage,
      perFamilyGidToSpid,
      distanceMatrices);

  unsigned int speciesNumber = mappings.getCoveredSpecies().size();
  Logger::info << "Species: " << speciesNumber << std::endl;
  Logger::info << "Distance matrices: " << distanceMatrices.size() << std::endl;

  auto speciesLabels = mappings.getCoveredSpecies();
  SpeciesTrees startingSpeciesTrees =
    generateRandomSpeciesTrees(speciesLabels, 
        speciesToSpid,
        1);
  auto &speciesTree = *startingSpeciesTrees[0];


  AsteroidOptimizer optimizer(speciesTree,
      perFamilyCoverage,
      perFamilyGidToSpid,
      distanceMatrices);
  optimizer.optimize();
  std::string outputSpeciesTreeFile = arg.prefix + ".bestTree.newick";
  speciesTree.save(outputSpeciesTreeFile);
  close();
  return 0;

}
