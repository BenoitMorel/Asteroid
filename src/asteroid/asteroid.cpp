#include <iostream>
#include <corax/corax.h>

#include "AsteroidOptimizer.hpp"
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <maths/Random.hpp>
#include "Arguments.hpp"
#include <trees/PLLUnrootedTree.hpp>
#include <DistanceMethods/InternodeDistance.hpp>
#include <memory>
#include <set>

using GeneTrees = std::vector<std::shared_ptr<PLLUnrootedTree> >;

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
    const GeneSpeciesMapping &mapping,
    const StringToUint &speciesToSpid)
{
  for (auto leaf: geneTree.getLeaves()) {
    const auto &species = mapping.getSpecies(leaf->label);
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

int main(int argc, char * argv[])
{
  Logger::init();
  Logger::timed << "Starting Asteroid..." << std::endl;
  Arguments arg(argc, argv);
  arg.printCommand();
  Random::setSeed(static_cast<unsigned int>(arg.seed));

  // read all gene trees
  Logger::timed << "Parsing gene trees..." << std::endl;
  GeneTrees geneTrees;
  parseTrees(arg.inputGeneTreeFile,
      geneTrees);

  Logger::timed << "Number of gene trees:" << geneTrees.size() << std::endl;

  // fill gene to species mapping
  GeneSpeciesMapping mapping;
  if (!arg.geneSpeciesMapping.size()) {
    for (auto geneTree: geneTrees) {
      mapping.fillFromGeneTree(*geneTree);
    }
  } else {
    mapping.fillFromMappingFile(arg.geneSpeciesMapping); 
  }
 
  // species to ID mapping
  StringToUint speciesToSpid;
  unsigned int speciesNumber = mapping.getCoveredSpecies().size();
  for (const auto &species: mapping.getCoveredSpecies()) {
    speciesToSpid.insert({species, speciesToSpid.size()});
  }

  Logger::timed << "Number of species: " 
    << speciesNumber
    << std::endl;

  for (auto geneTree: geneTrees) {
    setGeneTreeSpid(*geneTree, mapping, speciesToSpid);
  }

  Logger::timed << "Computing internode distances from the gene trees..." << std::endl;
  std::vector<DistanceMatrix> distanceMatrices;
  IDParam idParams;
  idParams.minBL = arg.minBL;
 
  if (!arg.noCorrection) {
    fillDistanceMatricesAsteroid(geneTrees, 
        speciesNumber, 
        idParams, 
        distanceMatrices);
  } else {
    fillDistanceMatricesAstrid(geneTrees, 
        speciesNumber, 
        idParams, 
        distanceMatrices);
  }
  Logger::timed << "Initializing random starting species tree..." << std::endl;
  // Initial species tree 
  PLLUnrootedTree speciesTree(mapping.getCoveredSpecies());
  speciesTree.reindexLeaves(speciesToSpid);
  for (auto leaf: speciesTree.getLeaves()) {
    intptr_t spid = static_cast<intptr_t>(leaf->node_index);
    leaf->data = reinterpret_cast<void*>(spid);
  }

    

  BoolMatrix perFamilyCoverage;
  
  if (!arg.noCorrection) {
    perFamilyCoverage = BoolMatrix(geneTrees.size(),
        std::vector<bool>(speciesNumber, false));
    for (unsigned int k = 0; k < geneTrees.size(); ++k) {
      for (auto geneLeaf: geneTrees[k]->getLeaves()) {
        auto spid = reinterpret_cast<intptr_t>(geneLeaf->data);
        perFamilyCoverage[k][spid] = true;
      }
    }
  } else {
    perFamilyCoverage.push_back(std::vector<bool>(speciesNumber, true));
  }
  

  Logger::timed << "Initializing optimizer... " << std::endl;
  // Perform search
  AsteroidOptimizer optimizer(speciesTree,
      perFamilyCoverage,
      distanceMatrices);
  Logger::timed << "Starting tree search... " << std::endl;
  optimizer.optimize();

  std::string outputSpeciesTreeFile = arg.prefix + ".newick";
  speciesTree.save(outputSpeciesTreeFile);

  return 0;

}
