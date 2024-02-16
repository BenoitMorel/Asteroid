#include "GeneSpeciesMapping.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <IO/Logger.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <corax/corax.h>

void GeneSpeciesMapping::fill(const GeneSpeciesMapping &mapping)
{
  const auto &m = mapping.getMap();
  for (const auto &pair: m) {
    _map[pair.first] = pair.first;
  }
}
  

void GeneSpeciesMapping::fillFromMappingFile(const std::string &mappingFile)
{
  std::ifstream f(mappingFile);
  std::string line;
  getline(f, line);
  f.seekg(0, std::ios::beg);
  if (line.find(":") != std::string::npos) {
    buildFromPhyldogMapping(f);
  } else {
    buildFromTreerecsMapping(f);
  }
}

static bool isBlanck(const std::string &s) 
{
  return s.empty() 
    || std::all_of(s.begin(), 
        s.end(), 
        [](char c){return std::isspace(c);});
}

  

void GeneSpeciesMapping::buildFromPhyldogMapping(std::ifstream &f)
{
  /*
    species1:gene1;gene2;gene3
    species2:gene4;gene5
  */
  std::string line;
  while (getline(f, line)) {
    if (isBlanck(line)) {
      continue;
    }
    std::stringstream ss(line);
    std::string species;
    std::string gene;
    getline(ss, species, ':');
    _species.insert(species);
    while(getline(ss, gene, ';')) {
      _map[gene] = species;
    }
  }
}

void GeneSpeciesMapping::buildFromTreerecsMapping(std::ifstream &f)
{
  /*
    gene1 species1
    gene2 species2
    gene3 species1
  */
  std::string line;
  while (getline(f, line)) {
    if (isBlanck(line)) {
      continue;
    }
    std::stringstream ss(line);
    std::string species;
    std::string gene;
    ss >> gene;
    ss >> species;
    _map[gene] = species;
    _species.insert(species);
  }
}


void GeneSpeciesMapping::fillFromGeneTree(PLLUnrootedTree &geneTree, bool fixDuplicates)
{
  unsigned int i = 0;
  for (auto leaf: geneTree.getLeaves()) {
    std::string geneLabel(leaf->label);
    auto speciesLabel = geneLabel;
    if (fixDuplicates) {
      while (geneExists(geneLabel)) {
        // duplicated gene name, we rename it
        geneLabel = speciesLabel + "_" + std::to_string(i);
        ++i;
      }
      if (geneLabel != speciesLabel) {
        // the gene name has chaned, update it
        geneTree.setLabel(leaf->node_index, geneLabel);
      }
    } else {
      if (geneExists(geneLabel)) {
        std::cerr << "Warning, duplicated gene " + geneLabel + " in one of the gene trees. Aborting." << std::endl;
        assert(false);
      }
    }
    _map[geneLabel] = std::string(speciesLabel);
    _species.insert(speciesLabel);
  }
}

