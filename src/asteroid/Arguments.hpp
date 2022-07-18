#pragma once

#include <string>

struct Arguments {
  Arguments(int argc, char * argv[]);
  void printHelp();
  void printCommand();
 
  int argc;
  char ** argv;

  std::string inputGeneTreeFile;
  std::string inputBSGeneTreeFile;
  std::string inputWeights;
  unsigned int randomStartingTrees;
  unsigned int bootstrapReplicates;
  std::string prefix;
  
  // if this argument is not set, we assume
  // single-copy families and that gene names
  // correspond to species names
  std::string geneSpeciesMapping;
  unsigned int seed;

  bool noCorrection;
  double minBL;
  bool useGeneBL;
  bool stepwise;
};
