#pragma once

#include <string>
#include <map>
#include <fstream>
#include <unordered_set>
#include <unordered_map>

typedef struct corax_utree_s corax_utree_t;
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;
class PLLUnrootedTree;


/**
 *  Parse and stores gene-to-species mapping
 *
 */
class GeneSpeciesMapping {
public:
  
  GeneSpeciesMapping() {}

  /**
   *  Add all labels from the gene tree
   *  Each gene label should be equal to its species name
   */
  void fillFromGeneTree(PLLUnrootedTree &geneTree); 
  void fill(const GeneSpeciesMapping &mapping);
 
  /*
   *  Support the following format:
   *    species1:gene1;gene2;gene3
   *    species2:gene1
   *    species3:gene1:gene2;gene3
   *
   *    gene1 species1
   *    gene2 species2
   *    gene3 species1
  */
  void fillFromMappingFile(const std::string &mappingFile);

  /**
   *  @return an object mapping each gene to one species
   */
  const std::unordered_map<std::string, std::string> &getMap() const {return _map;}

  /**
   *  @param gene gene std::string
   *  @return the species mapped to this gene
   */
  const std::string &getSpecies(const std::string &gene) const {
    return _map.find(gene)->second;
  }

  const std::unordered_set<std::string> &getCoveredSpecies() const
  {
    return _species;
  }
private:
  std::unordered_map<std::string, std::string> _map; // <gene,species>
  std::unordered_set<std::string> _species;
  void buildFromPhyldogMapping(std::ifstream &f);
  void buildFromTreerecsMapping(std::ifstream &f);
};

